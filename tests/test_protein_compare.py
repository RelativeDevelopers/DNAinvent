"""Tests de la fase 4 de proteínas: superposición, comparación, ESMFold y predicción en el
recorrido completo. Sin red: ESMFold recibe una descarga falsa."""

import math
import random
import urllib.error
from dataclasses import replace
from pathlib import Path

import pytest

from dnainvent import InvalidParameterError, ServiceError
from dnainvent.protein import confidence, pipeline
from dnainvent.protein.cache import Cache
from dnainvent.protein.classify import Category
from dnainvent.protein.clients.blast import SWISSPROT, BlastHit
from dnainvent.protein.compare import (Difference, compare, moved, pairs_by_alignment,
                                       pairs_from_blast)
from dnainvent.protein.pdb import Structure, read_pdb, write_pdb
from dnainvent.protein.predictors.esmfold import API_URL, ESMFold
from dnainvent.protein.superpose import rmsd, superpose
from test_protein_search import HBB, FakeAlphaFold, FakeBlast, blast_hit, gene

DATA = Path(__file__).parent / "data"
SICKLE = HBB.replace("MVHLTPEEK", "MVHLTPVEK", 1)        # Glu7Val (numeración de UniProt)


@pytest.fixture(scope="module")
def hemoglobin():
    return read_pdb(DATA / "AF-P68871-F1-model_v6.pdb")


def random_rotation(rng):
    q = [rng.gauss(0, 1) for _ in range(4)]
    n = math.sqrt(sum(c * c for c in q))
    w, x, y, z = (c / n for c in q)
    return ((w*w + x*x - y*y - z*z, 2*(x*y - w*z), 2*(x*z + w*y)),
            (2*(x*y + w*z), w*w - x*x + y*y - z*z, 2*(y*z - w*x)),
            (2*(x*z - w*y), 2*(y*z + w*x), w*w - x*x - y*y + z*z))


def rotate(r, t, p):
    return tuple(sum(r[i][k] * p[k] for k in range(3)) + t[i] for i in range(3))


def as_prediction(structure: Structure, sequence: str, rng=None, plddt=None) -> Structure:
    """La estructura, girada y desplazada, con los aminoácidos de 'sequence' (lo que
    devolvería un predictor que acierta la forma)."""
    from dnainvent.constants import AMINO_ACID_CODES
    rng = rng or random.Random(1)
    r, t = random_rotation(rng), tuple(rng.uniform(-30, 30) for _ in range(3))
    index = {resi: k for k, resi in enumerate(sorted({a.resi for a in structure.atoms}))}
    atoms = []
    for a in structure.atoms:
        x, y, z = rotate(r, t, a.xyz)
        resn = AMINO_ACID_CODES[sequence[index[a.resi]]].upper()
        atoms.append(replace(a, x=x, y=y, z=z, resn=resn,
                             bfactor=a.bfactor if plddt is None else plddt))
    return Structure(atoms, "ESMFOLD V1 PREDICTION FOR INPUT")


# ---------- Superposición ----------

def test_superpose_recovers_a_rigid_motion():
    rng = random.Random(4)
    for _ in range(50):
        points = [tuple(rng.uniform(-20, 20) for _ in range(3)) for _ in range(rng.randint(3, 40))]
        r, t = random_rotation(rng), tuple(rng.uniform(-50, 50) for _ in range(3))
        target = [rotate(r, t, p) for p in points]
        fit = superpose(points, target)
        assert fit.rmsd < 1e-9
        assert rmsd([fit.apply(p) for p in points], target) < 1e-9


def test_superpose_with_noise_gives_the_noise_level():
    rng = random.Random(8)
    points = [tuple(rng.uniform(-20, 20) for _ in range(3)) for _ in range(500)]
    noisy = [tuple(c + rng.gauss(0, 0.5) for c in p) for p in points]
    assert 0.75 < superpose(points, noisy).rmsd < 0.95       # ~0,5·√3 ≈ 0,87 Å


def test_superpose_needs_three_pairs():
    with pytest.raises(ValueError):
        superpose([(0, 0, 0), (1, 0, 0)], [(0, 0, 0), (1, 0, 0)])


def test_superpose_never_mirrors():
    """Una imagen especular no se puede conseguir girando: el RMSD no puede ser 0."""
    rng = random.Random(2)
    points = [tuple(rng.uniform(-10, 10) for _ in range(3)) for _ in range(30)]
    mirror = [(-x, y, z) for x, y, z in points]
    fit = superpose(points, mirror)
    assert fit.rmsd > 1.0
    r = fit.rotation
    det = (r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
           - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
           + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]))
    assert det == pytest.approx(1.0)


# ---------- Comparación ----------

def test_sickle_cell_mutation(hemoglobin):
    c = compare(SICKLE, HBB, natural_structure=hemoglobin,
                predicted=as_prediction(hemoglobin, SICKLE))
    assert [d.label for d in c.differences] == ["Glu7Val"]
    assert c.identity == pytest.approx(146 / 147 * 100)
    assert c.rmsd < 1e-6 and c.superposition.count == 147
    assert max(c.deviations.values()) < 1e-6
    assert len(c.conserved) == 146


def test_moved_prediction_lands_on_the_natural(hemoglobin):
    predicted = as_prediction(hemoglobin, HBB, random.Random(9))
    c = compare(HBB, HBB, natural_structure=hemoglobin, predicted=predicted)
    placed = moved(predicted, c.superposition)
    assert max(math.dist(a.xyz, b.xyz) for a, b in zip(placed.atoms, hemoglobin.atoms)) < 1e-6


def test_pairs_from_blast_with_gaps():
    user = HBB[:20] + HBB[25:]                              # le faltan 5 aminoácidos
    hit = replace(blast_hit(), query_from=1, query_to=len(user), hit_from=1, hit_to=147,
                  query_aligned=HBB[:20] + "-----" + HBB[25:], hit_aligned=HBB)
    pairs = pairs_from_blast(hit)
    assert len(pairs) == 142 and pairs[19:21] == [(19, 19), (20, 25)]
    c = compare(user, HBB, hit=hit)
    assert (c.insertions, c.deletions, c.differences) == (0, 5, [])


def test_pairs_from_blast_start_in_the_middle():
    hit = replace(blast_hit(), query_from=3, hit_from=11, query_aligned="AB-C", hit_aligned="AXDC")
    assert pairs_from_blast(hit) == [(2, 10), (3, 11), (4, 13)]


def test_without_blast_alignment_uses_our_own():
    user = "MKT" + HBB[3:]
    c = compare(user, HBB)
    assert [d.label for d in c.differences] == ["Val2Lys", "His3Thr"]
    assert pairs_by_alignment(user, HBB)[:3] == [(0, 0), (1, 1), (2, 2)]


def test_difference_label_uses_three_letter_codes():
    assert Difference(7, 7, "E", "V").label == "Glu7Val"
    assert Difference(10, 9, "W", "X").label == "Trp10Xaa"


# ---------- ESMFold (sin red) ----------

def esmfold_answer(structure: Structure, scale: float = 1.0) -> bytes:
    atoms = [replace(a, bfactor=round(a.bfactor * scale, 4)) for a in structure.atoms]
    return write_pdb(atoms, "ESMFOLD V1 PREDICTION FOR INPUT").encode()


class FakeESMAtlas:
    def __init__(self, answer):
        self.answer, self.calls = answer, []

    def __call__(self, url, text):
        self.calls.append((url, text))
        if isinstance(self.answer, Exception):
            raise self.answer
        return self.answer


def test_esmfold_prediction(tmp_path, hemoglobin):
    post = FakeESMAtlas(esmfold_answer(hemoglobin))
    structure = ESMFold(Cache(tmp_path), post).predict(HBB.lower())
    assert post.calls == [(API_URL, HBB)]
    assert structure.sequence() == HBB and structure.is_prediction
    assert confidence.mean(structure) == pytest.approx(confidence.mean(hemoglobin))


def test_esmfold_plddt_from_0_to_1_is_scaled(tmp_path, hemoglobin):
    post = FakeESMAtlas(esmfold_answer(hemoglobin, scale=0.01))
    structure = ESMFold(Cache(tmp_path), post).predict(HBB)
    # el PDB guarda 2 decimales: 0,9718 se escribe 0,97, de ahí la tolerancia
    assert confidence.mean(structure) == pytest.approx(confidence.mean(hemoglobin), abs=0.5)
    assert max(a.bfactor for a in structure.atoms) > 1


def test_esmfold_uses_the_cache(tmp_path, hemoglobin):
    ESMFold(Cache(tmp_path), FakeESMAtlas(esmfold_answer(hemoglobin))).predict(HBB)
    again = FakeESMAtlas(b"no se debe pedir")
    assert ESMFold(Cache(tmp_path), again).predict(HBB).sequence() == HBB
    assert again.calls == []


def test_esmfold_too_long(tmp_path):
    post = FakeESMAtlas(b"")
    with pytest.raises(InvalidParameterError) as info:
        ESMFold(Cache(tmp_path), post).predict("M" * 401)
    assert info.value.code == "prediction_too_long"
    assert info.value.params == {"length": 401, "maximum": 400}
    assert post.calls == []


@pytest.mark.parametrize("answer, code", [
    (urllib.error.URLError("certificado SSL"), "service_unavailable"),
    (urllib.error.HTTPError(API_URL, 503, "x", None, None), "service_unavailable"),
    (TimeoutError(), "service_unavailable"),
    (urllib.error.HTTPError(API_URL, 400, "x", None, None), "prediction_failed"),  # rechazada
    (b"<html>Internal Server Error</html>", "prediction_failed"),
    (b"", "prediction_failed"),
])
def test_esmfold_failures(tmp_path, answer, code):
    with pytest.raises(ServiceError) as info:
        ESMFold(Cache(tmp_path), FakeESMAtlas(answer)).predict(HBB)
    assert info.value.code == code


def test_esmfold_wrong_sequence_is_not_cached(tmp_path, hemoglobin):
    with pytest.raises(ServiceError):
        ESMFold(Cache(tmp_path), FakeESMAtlas(esmfold_answer(hemoglobin))).predict(SICKLE)
    assert not list(tmp_path.rglob("*.pdb"))


# ---------- Recorrido completo con predicción ----------

class FakePredictor:
    name, max_length = "ESMFold", 400

    def __init__(self, structure=None, error=None, plddt=None):
        self.structure, self.error, self.plddt, self.calls = structure, error, plddt, []

    def predict(self, protein):
        self.calls.append(protein)
        if self.error:
            raise self.error
        return as_prediction(self.structure, protein, plddt=self.plddt)


def test_identical_protein_is_not_predicted(hemoglobin):
    predictor = FakePredictor(hemoglobin)
    report = pipeline.analyze(gene(HBB), FakeBlast({SWISSPROT: [blast_hit()]}), FakeAlphaFold(),
                              predictor)
    assert predictor.calls == [] and report.prediction is None and report.comparison is None


def test_variant_is_predicted_and_compared(hemoglobin):
    hit = replace(blast_hit(identity=99.3), query_aligned=SICKLE, hit_aligned=HBB,
                  markup=HBB.replace("E", " ", 1))
    predictor, stages = FakePredictor(hemoglobin), []
    report = pipeline.analyze(gene(SICKLE), FakeBlast({SWISSPROT: [hit]}), FakeAlphaFold(),
                              predictor, progress=lambda key, data: stages.append(key))
    assert report.category is Category.VERY_SIMILAR
    assert predictor.calls == [SICKLE] and report.predictor == "ESMFold"
    assert [d.label for d in report.comparison.differences] == ["Glu7Val"]
    assert report.comparison.rmsd < 1e-6
    assert stages[-1] == "protein.stage.predict"


def test_variant_without_natural_structure_still_lists_mutations(hemoglobin):
    hit = replace(blast_hit(identity=99.3), query_aligned=SICKLE, hit_aligned=HBB)
    afdb = FakeAlphaFold(ServiceError("protein_not_found", accession="P68871"))
    report = pipeline.analyze(gene(SICKLE), FakeBlast({SWISSPROT: [hit]}), afdb,
                              FakePredictor(hemoglobin))
    assert [d.label for d in report.comparison.differences] == ["Glu7Val"]
    assert report.comparison.rmsd is None


def test_invented_protein_is_predicted(hemoglobin):
    predictor = FakePredictor(hemoglobin, plddt=35.0)
    report = pipeline.analyze(gene(HBB), FakeBlast({}), FakeAlphaFold(), predictor)
    assert report.category is Category.INVENTED and predictor.calls == [HBB]
    assert report.comparison is None
    assert ("protein.warnings.unlikely_fold", {"plddt": 35.0}) in report.warnings


def test_invented_protein_with_good_confidence_has_no_fold_warning(hemoglobin):
    report = pipeline.analyze(gene(HBB), FakeBlast({}), FakeAlphaFold(),
                              FakePredictor(hemoglobin, plddt=85.0))
    assert not any(key == "protein.warnings.unlikely_fold" for key, _ in report.warnings)


def test_short_proteins_are_not_predicted(hemoglobin):
    predictor = FakePredictor(hemoglobin)
    report = pipeline.analyze(gene(HBB[:12]), FakeBlast({}), None, predictor)
    assert predictor.calls == [] and report.prediction is None
    assert [key for key, _ in report.warnings] == ["protein.warnings.too_short"]


@pytest.mark.parametrize("error, warning", [
    (InvalidParameterError("prediction_too_long", length=500, maximum=400),
     ("protein.warnings.prediction_too_long", {"length": 500, "maximum": 400})),
    (ServiceError("service_unavailable", service="ESMFold"),
     ("protein.warnings.prediction_unavailable", {"predictor": "ESMFold"})),
    (ServiceError("prediction_failed", service="ESMFold"),
     ("protein.warnings.prediction_failed", {"predictor": "ESMFold"})),
])
def test_prediction_problems_do_not_stop_the_report(error, warning):
    report = pipeline.analyze(gene(HBB), FakeBlast({}), None, FakePredictor(error=error))
    assert report.category is Category.INVENTED and report.prediction is None
    assert warning in report.warnings
