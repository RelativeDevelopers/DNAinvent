"""Tests de la fase 3 de proteínas: cliente BLAST del EBI y recorrido completo (pipeline).

Sin red: BLAST y AlphaFold DB reciben descargas falsas.
"""

import json
import urllib.error
from pathlib import Path

import pytest

from dnainvent import InvalidParameterError, InvalidSequenceError, ServiceError
from dnainvent.protein import pipeline
from dnainvent.protein.cache import Cache
from dnainvent.protein.classify import Category
from dnainvent.protein.clients.afdb import AlphaFoldEntry
from dnainvent.protein.clients import blast as blast_module
from dnainvent.protein.clients.blast import (SWISSPROT, UNIPROTKB, BlastClient, BlastHit,
                                             BlastResult, parse_result, validate_email)
from dnainvent.protein.orf import best_orf
from dnainvent.protein.pdb import read_pdb

DATA = Path(__file__).parent / "data"
HBB = ("MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDN"
       "LKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH")
EMAIL = "persona@example.org"
BASE = blast_module.BASE_URL

# Un codón para cada aminoácido: con esto se "retrotraduce" una proteína a ADN.
CODON = {"A": "GCT", "R": "CGT", "N": "AAC", "D": "GAC", "C": "TGC", "Q": "CAG", "E": "GAA",
         "G": "GGT", "H": "CAC", "I": "ATC", "L": "CTG", "K": "AAA", "M": "ATG", "F": "TTC",
         "P": "CCG", "S": "TCT", "T": "ACC", "W": "TGG", "Y": "TAC", "V": "GTT"}


def gene(protein: str, stop: str = "TAA") -> str:
    return "".join(CODON[aa] for aa in protein) + stop


def ebi_hit(accession="P68871", identity=100.0, q=(1, 147), h=(1, 147), length=147,
            evalue=1e-100, desc="Hemoglobin subunit beta OS=Homo sapiens OX=9606 GN=HBB PE=1 SV=2",
            entry="HBB_HUMAN", uni=True):
    """Un hit con el formato JSON del EBI."""
    item = {
        "hit_num": 1, "hit_db": "SP", "hit_id": entry, "hit_acc": accession,
        "hit_desc": desc, "hit_len": length,
        "hit_hsps": [{
            "hsp_num": 1, "hsp_score": 760, "hsp_bit_score": 297.4, "hsp_expect": evalue,
            "hsp_align_len": q[1] - q[0] + 1, "hsp_identity": identity, "hsp_positive": 100.0,
            "hsp_gaps": 0, "hsp_query_from": q[0], "hsp_query_to": q[1],
            "hsp_hit_from": h[0], "hsp_hit_to": h[1],
            "hsp_qseq": "MVHLT", "hsp_mseq": "MVHLT", "hsp_hseq": "MVHLT",
        }],
    }
    if uni:
        item.update({"hit_uni_de": "Hemoglobin subunit beta", "hit_uni_os": "Homo sapiens",
                     "hit_uni_gn": "HBB"})
    return item


def ebi_json(*hits) -> bytes:
    return json.dumps({"program": "blastp", "query_len": 147, "hits": list(hits)}).encode()


class FakeEBI:
    """Simula el Job Dispatcher: /run, /status (con una lista de estados) y /result."""

    def __init__(self, result: bytes = b"", statuses=("RUNNING", "FINISHED"), run_error=None):
        self.result, self.statuses = result, list(statuses)
        self.run_error = run_error
        self.posts, self.gets, self.slept = [], [], []
        self.time = 0.0

    def post(self, url, data):
        self.posts.append((url, data))
        if self.run_error:
            raise self.run_error
        return b"ncbiblast-R20261102-000001-0001-1234-p1m"

    def get(self, url):
        self.gets.append(url)
        if "/status/" in url:
            status = self.statuses.pop(0) if len(self.statuses) > 1 else self.statuses[0]
            if isinstance(status, Exception):
                raise status
            return status.encode()
        return self.result

    def sleep(self, seconds):
        self.slept.append(seconds)
        self.time += seconds

    def clock(self):
        return self.time

    def client(self, tmp_path, email=EMAIL, timeout=blast_module.TIMEOUT):
        return BlastClient(email, Cache(tmp_path), self.get, self.post, self.sleep, self.clock,
                           timeout=timeout)


# ---------- Cliente BLAST ----------

def test_blast_full_job(tmp_path):
    ebi = FakeEBI(ebi_json(ebi_hit()), statuses=["QUEUED", "RUNNING", "FINISHED"])
    seen = []
    result = ebi.client(tmp_path).search(HBB, progress=seen.append)
    url, data = ebi.posts[0]
    assert url == f"{BASE}/run"
    assert data == {"email": EMAIL, "program": "blastp", "stype": "protein",
                    "sequence": HBB, "database": SWISSPROT}
    assert seen == ["QUEUED", "RUNNING"]
    assert ebi.gets[-1].endswith("/json")
    assert ebi.slept[0] == blast_module.FIRST_WAIT
    (hit,) = result.hits
    assert (hit.accession, hit.entry_name, hit.description, hit.organism, hit.gene) == \
        ("P68871", "HBB_HUMAN", "Hemoglobin subunit beta", "Homo sapiens", "HBB")
    assert (hit.identity, hit.query_coverage, hit.subject_coverage) == (100, 100, 100)


def test_blast_results_are_cached(tmp_path):
    ebi = FakeEBI(ebi_json(ebi_hit()))
    ebi.client(tmp_path).search(HBB)
    again = FakeEBI(b"no se debe pedir")
    assert again.client(tmp_path).search(HBB).hits[0].accession == "P68871"
    assert again.posts == [] and again.gets == []
    other_db = FakeEBI(ebi_json())
    assert other_db.client(tmp_path).search(HBB, UNIPROTKB).hits == []    # otra base, otra clave
    assert other_db.posts


def test_description_parsed_from_uniprot_header():
    raw = ebi_json(ebi_hit(uni=False, desc="Myoglobin OS=Physeter catodon OX=9755 GN=MB PE=1 SV=2"))
    hit = parse_result(raw, HBB, SWISSPROT).hits[0]
    assert (hit.description, hit.organism, hit.gene) == ("Myoglobin", "Physeter catodon", "MB")


def test_coverage_is_computed_from_the_aligned_stretch():
    raw = ebi_json(ebi_hit(q=(11, 110), h=(51, 150), length=400, identity=45.5))
    hit = parse_result(raw, "M" * 200, SWISSPROT).hits[0]
    assert hit.query_coverage == pytest.approx(50.0)       # 100 de los 200 de la consulta
    assert hit.subject_coverage == pytest.approx(25.0)     # 100 de los 400 de la real
    assert hit.to_hit().identity == 45.5


def test_hits_without_alignments_are_ignored():
    raw = json.dumps({"hits": [{"hit_acc": "X", "hit_hsps": []}, ebi_hit()]})
    assert [h.accession for h in parse_result(raw, HBB, SWISSPROT).hits] == ["P68871"]


@pytest.mark.parametrize("statuses, code", [
    (["ERROR"], "search_failed"),
    (["FAILURE"], "search_failed"),
    (["NOT_FOUND"], "search_failed"),
    ([urllib.error.URLError("sin red")], "service_unavailable"),
    ([urllib.error.HTTPError("u", 503, "x", None, None)], "service_unavailable"),
])
def test_blast_job_errors(tmp_path, statuses, code):
    with pytest.raises(ServiceError) as info:
        FakeEBI(statuses=statuses).client(tmp_path).search(HBB)
    assert info.value.code == code


def test_blast_rejected_submission(tmp_path):
    ebi = FakeEBI(run_error=urllib.error.HTTPError("u", 400, "Bad Request", None, None))
    with pytest.raises(ServiceError) as info:
        ebi.client(tmp_path).search(HBB)
    assert info.value.code == "search_failed"


def test_blast_timeout(tmp_path):
    ebi = FakeEBI(statuses=["RUNNING"])
    with pytest.raises(ServiceError) as info:
        ebi.client(tmp_path, timeout=60).search(HBB)
    assert info.value.code == "search_timeout" and info.value.params["minutes"] == 1


def test_blast_garbage_result(tmp_path):
    with pytest.raises(ServiceError) as info:
        FakeEBI(b"<html>error</html>").client(tmp_path).search(HBB)
    assert info.value.code == "search_failed"


@pytest.mark.parametrize("email", [None, "", "sin-arroba", "a@b", "dos@@x.com"])
def test_email_is_required_before_any_request(tmp_path, email):
    ebi = FakeEBI(ebi_json())
    with pytest.raises(InvalidParameterError) as info:
        ebi.client(tmp_path, email=email).search(HBB)
    assert info.value.code == "invalid_email"
    assert ebi.posts == []


def test_valid_email():
    assert validate_email("  persona@example.org ") == "persona@example.org"


def test_unknown_database(tmp_path):
    with pytest.raises(ValueError):
        FakeEBI().client(tmp_path).search(HBB, "pdb")


# ---------- Recorrido completo ----------

def blast_hit(accession="P68871", identity=100.0, query_coverage=100.0, subject_coverage=100.0,
              evalue=1e-100, length=147, hit_from=1, hit_to=147):
    return BlastHit(accession, f"{accession}_HUMAN", "Hemoglobin subunit beta", "Homo sapiens",
                    "HBB", length, identity, query_coverage, subject_coverage, evalue, 297.4,
                    1, 147, hit_from, hit_to, "MVHLT", "MVHLT", "MVHLT")


class FakeBlast:
    def __init__(self, answers):
        self.answers, self.calls = answers, []

    def search(self, protein, database=SWISSPROT, progress=None):
        self.calls.append((protein, database))
        if progress:
            progress("RUNNING")
        answer = self.answers.get(database, [])
        if isinstance(answer, Exception):
            raise answer
        return BlastResult(database, protein, answer)


class FakeAlphaFold:
    def __init__(self, error=None):
        self.error, self.calls = error, []

    def structure(self, accession):
        self.calls.append(accession)
        if self.error:
            raise self.error
        entry = AlphaFoldEntry(accession, f"AF-{accession}-F1", 6, HBB, "https://x/y.pdb", 97.19)
        return entry, read_pdb(DATA / "AF-P68871-F1-model_v6.pdb")


def test_the_back_translated_gene_gives_back_the_protein():
    assert best_orf("CC" + gene(HBB) + "G").protein == HBB


def test_known_protein(tmp_path):
    blast, afdb, stages = FakeBlast({SWISSPROT: [blast_hit()]}), FakeAlphaFold(), []
    report = pipeline.analyze(gene(HBB), blast, afdb,
                              progress=lambda key, data: stages.append(key))
    assert report.category is Category.IDENTICAL
    assert report.protein == HBB and report.best.accession == "P68871"
    assert report.database == SWISSPROT and blast.calls == [(HBB, SWISSPROT)]
    assert afdb.calls == ["P68871"] and report.structure.sequence() == HBB
    assert report.warnings == []
    assert stages == ["protein.stage.orf", "protein.stage.search", "protein.stage.waiting",
                      "protein.stage.structure"]


def test_variant_warns_about_mutations():
    report = pipeline.analyze(gene(HBB), FakeBlast({SWISSPROT: [blast_hit(identity=99.3)]}),
                              FakeAlphaFold())
    assert report.category is Category.VERY_SIMILAR
    assert ("protein.warnings.mutation_effect", {}) in report.warnings


def test_search_continues_in_uniprotkb():
    blast = FakeBlast({SWISSPROT: [], UNIPROTKB: [blast_hit("A0A0A0", identity=60)]})
    report = pipeline.analyze(gene(HBB), blast, FakeAlphaFold())
    assert [db for _, db in blast.calls] == [SWISSPROT, UNIPROTKB]
    assert report.category is Category.DISTANT_RELATIVE and report.database == UNIPROTKB


def test_invented_protein():
    blast, afdb = FakeBlast({}), FakeAlphaFold()
    report = pipeline.analyze(gene(HBB), blast, afdb)
    assert report.category is Category.INVENTED and report.best is None
    assert [db for _, db in blast.calls] == [SWISSPROT, UNIPROTKB]
    assert afdb.calls == []
    assert report.prediction is None                       # sin predictor no se predice


def test_only_swissprot_if_asked():
    blast = FakeBlast({})
    pipeline.analyze(gene(HBB), blast, None, search_everything=False)
    assert [db for _, db in blast.calls] == [SWISSPROT]


def test_unreliable_hits_do_not_count():
    weak = blast_hit(identity=100, query_coverage=30)           # cubre poco de la consulta
    report = pipeline.analyze(gene(HBB), FakeBlast({SWISSPROT: [weak], UNIPROTKB: [weak]}))
    assert report.category is Category.INVENTED
    assert report.hits == [weak]                                 # se enseñan igualmente


def test_too_short_to_search():
    blast = FakeBlast({SWISSPROT: [blast_hit()]})
    report = pipeline.analyze(gene("MKVLAT"), blast)
    assert blast.calls == [] and report.category is Category.INVENTED
    assert ("protein.warnings.too_short", {"length": 6, "minimum": 20}) in report.warnings


def test_incomplete_reading_frame_is_reported():
    report = pipeline.analyze(gene(HBB, stop=""), FakeBlast({SWISSPROT: [blast_hit()]}),
                              FakeAlphaFold())
    assert ("protein.warnings.incomplete_orf", {}) in report.warnings


def test_fragment_of_a_real_protein():
    fragment = blast_hit(identity=100, subject_coverage=40, length=370, hit_from=21, hit_to=167)
    report = pipeline.analyze(gene(HBB), FakeBlast({SWISSPROT: [fragment]}), FakeAlphaFold())
    assert report.category is Category.VERY_SIMILAR
    assert ("protein.warnings.fragment", {"start": 21, "end": 167, "length": 370}) in report.warnings


@pytest.mark.parametrize("code, warning", [
    ("protein_not_found", "protein.warnings.no_structure"),
    ("service_unavailable", "protein.warnings.structure_unavailable"),
])
def test_structure_problems_do_not_stop_the_report(code, warning):
    params = {"accession": "P68871"} if code == "protein_not_found" else {"service": "AlphaFold DB"}
    afdb = FakeAlphaFold(ServiceError(code, **params))
    report = pipeline.analyze(gene(HBB), FakeBlast({SWISSPROT: [blast_hit()]}), afdb)
    assert report.category is Category.IDENTICAL and report.structure is None
    assert (warning, {"accession": "P68871"}) in report.warnings


def test_blast_failure_reaches_the_interface():
    blast = FakeBlast({SWISSPROT: ServiceError("search_failed", service="BLAST (EBI)")})
    with pytest.raises(ServiceError):
        pipeline.analyze(gene(HBB), blast)


def test_no_reading_frame():
    with pytest.raises(InvalidSequenceError) as info:
        pipeline.analyze("CCCCCCCCCC", FakeBlast({}))
    assert info.value.code == "no_orf"


def test_every_warning_and_stage_has_texts():
    """Cada clave protein.stage.* / protein.warnings.* que usa pipeline.py, y el nombre de
    cada base de datos, tiene texto en todos los idiomas."""
    import re
    from i18n import Translator, available_languages
    source = Path(pipeline.__file__).read_text(encoding="utf-8")
    keys = set(re.findall(r'"(protein\.(?:stage|warnings)\.[a-z_]+)"', source))
    keys |= {f"protein.databases.{db}" for db in blast_module.DATABASES}
    assert len(keys) >= 13
    for language in available_languages():
        tr = Translator(language)
        for key in keys:
            assert tr(key) != key, (language, key)
