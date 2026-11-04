"""Tests de la fase 2 de proteínas: lector PDB, hidrógenos, molécula 3D, AlphaFold DB y caché.

Sin red: el cliente de AlphaFold DB recibe una función de descarga falsa.
"""

import json
import os
import random
import re
import time
import urllib.error
from dataclasses import replace
from pathlib import Path
from statistics import mean, median

import pytest

from dnainvent import InvalidParameterError, ServiceError
from dnainvent.protein import confidence
from dnainvent.protein.cache import Cache, key_for
from dnainvent.protein.clients.afdb import AlphaFoldDB, validate_accession
from dnainvent.protein.geometry import angle, dihedral, distance, place
from dnainvent.protein.hydrogens import add_hydrogens, is_rotatable, templates
from dnainvent.protein.pdb import normalize_name, parse_pdb, read_pdb, write_pdb

DATA = Path(__file__).parent / "data"
HEMOGLOBIN = DATA / "AF-P68871-F1-model_v6.pdb"
UBIQUITIN = DATA / "1d3z_model1.pdb"
HBB = ("MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDN"
       "LKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH")
UBQ = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"


@pytest.fixture(scope="module")
def hemoglobin():
    return read_pdb(HEMOGLOBIN)


@pytest.fixture(scope="module")
def ubiquitin():
    return read_pdb(UBIQUITIN)


def atom_line(name, resn="ALA", resi=1, x=0.0, alt=" ", record="ATOM  ", element=None):
    element = element if element is not None else name[0]
    return (f"{record}{1:5d} {name:<4}{alt}{resn:>3} A{resi:4d}    {x:8.3f}{0:8.3f}{0:8.3f}"
            f"{1:6.2f}{50:6.2f}          {element:>2}")


# ---------- Lector PDB ----------

def test_read_alphafold_file(hemoglobin):
    assert hemoglobin.is_prediction
    assert "HEMOGLOBIN SUBUNIT BETA" in hemoglobin.title
    assert hemoglobin.chains() == ["A"]
    assert len(hemoglobin.atoms) == 1131
    assert hemoglobin.sequence() == HBB
    ca = hemoglobin.residues()[0].atoms["CA"]
    assert (ca.resn, ca.resi, ca.element) == ("MET", 1, "C")
    assert ca.bfactor == pytest.approx(61.09)        # en AlphaFold, la columna B es el pLDDT


def test_read_nmr_file_with_hydrogens(ubiquitin):
    assert not ubiquitin.is_prediction
    assert ubiquitin.sequence() == UBQ
    assert sum(a.element == "H" for a in ubiquitin.atoms) == 629


def test_only_the_requested_model():
    text = "\n".join(["MODEL        1", atom_line("CA", x=1.0), "ENDMDL",
                      "MODEL        2", atom_line("CA", x=2.0), "ENDMDL"])
    assert [a.x for a in parse_pdb(text).atoms] == [1.0]
    assert [a.x for a in parse_pdb(text, model=2).atoms] == [2.0]


def test_alternative_locations_keep_the_first():
    text = "\n".join([atom_line("CA", x=1.0, alt="A"), atom_line("CA", x=2.0, alt="B"),
                      atom_line("CB", x=3.0)])
    assert [(a.name, a.x) for a in parse_pdb(text).atoms] == [("CA", 1.0), ("CB", 3.0)]


def test_hetero_atoms_are_not_amino_acids():
    text = "\n".join([atom_line("CA"), atom_line("O", resn="HOH", resi=2, record="HETATM")])
    structure = parse_pdb(text)
    assert len(structure.atoms) == 2
    assert structure.sequence() == "A"


def test_element_from_name_when_column_is_missing():
    assert parse_pdb(atom_line("CA", element="")).atoms[0].element == "C"


def test_old_hydrogen_names_are_normalized():
    assert normalize_name("2HB") == "HB2"
    assert normalize_name("1HD2") == "HD21"
    assert normalize_name("HA") == "HA"


def test_write_and_read_back(hemoglobin):
    assert parse_pdb(write_pdb(hemoglobin.atoms)).atoms == hemoglobin.atoms


# ---------- Geometría ----------

def test_place_inverts_internal_coordinates():
    rng = random.Random(5)
    for _ in range(200):
        a, b, c, d = [tuple(rng.uniform(-5, 5) for _ in range(3)) for _ in range(4)]
        p = place(a, b, c, distance(c, d), angle(b, c, d), dihedral(a, b, c, d))
        assert distance(p, d) < 1e-9


# ---------- Hidrógenos ----------

# Hidrógenos de cada aminoácido dentro de la cadena, a pH 7.
H_PER_RESIDUE = {"ALA": 5, "ARG": 13, "ASN": 6, "ASP": 4, "CYS": 5, "GLN": 8, "GLU": 6,
                 "GLY": 3, "HIS": 7, "ILE": 11, "LEU": 11, "LYS": 13, "MET": 9, "PHE": 9,
                 "PRO": 7, "SER": 5, "THR": 7, "TRP": 10, "TYR": 9, "VAL": 9}


def test_hydrogen_counts_at_ph7():
    data = templates()["residues"]
    assert {resn: len(t["hydrogens"]) for resn, t in data.items()} == H_PER_RESIDUE


def test_protonation_choices():
    names = {resn: {h["name"] for h in t["hydrogens"]} for resn, t in templates()["residues"].items()}
    assert "HE2" in names["HIS"] and "HD1" not in names["HIS"]       # histidina con H en Nε
    assert {"HZ1", "HZ2", "HZ3"} <= names["LYS"]                      # Lys+
    assert not {"HD2"} & names["ASP"] and not {"HE2"} & names["GLU"]  # Asp-, Glu-
    assert "H" not in names["PRO"]                                    # la prolina no tiene H en N


def test_hemoglobin_gets_all_its_hydrogens(hemoglobin):
    result = add_hydrogens(hemoglobin.residues())
    expected = sum(H_PER_RESIDUE[r.resn] for r in hemoglobin.residues()) - 1 + 3  # extremo N
    assert result.added == expected == 1127
    assert result.skipped == []
    first = result.residues[0].atoms
    assert {"H1", "H2", "H3"} <= set(first) and "H" not in first      # -NH3+
    assert all(a.element == "H" for n, a in first.items() if n.startswith("H"))


def test_calculated_hydrogens_inherit_confidence(hemoglobin):
    residue = add_hydrogens(hemoglobin.residues()).residues[0]
    assert residue.atoms["HA"].bfactor == residue.atoms["CA"].bfactor


def compare_with_experiment(structure):
    residues = structure.residues()
    result = add_hydrogens(residues, keep_existing=False)
    fixed, methyl = [], []
    for original, calculated in zip(residues, result.residues):
        real = {n: a for n, a in original.atoms.items() if a.element == "H"}
        for name, atom in calculated.atoms.items():
            if atom.element != "H" or name not in real:
                continue
            if not is_rotatable(calculated.resn, name):
                fixed.append(distance(atom.xyz, real[name].xyz))
            elif re.fullmatch(r"H[BGDE]\d?[123]", name) and calculated.resn != "LYS":
                # en un metilo los tres H son intercambiables: el más cercano de los tres
                group = [a for n, a in real.items() if n[:-1] == name[:-1]]
                methyl.append(min(distance(atom.xyz, a.xyz) for a in group))
    return fixed, methyl


def test_fixed_hydrogens_match_the_real_ones(ubiquitin):
    """Ubiquitina 1D3Z (RMN): se quitan sus H, se recalculan y se comparan con los reales."""
    fixed, _ = compare_with_experiment(ubiquitin)
    assert len(fixed) > 400
    assert mean(fixed) < 0.07                          # ~0,05 Å
    assert max(fixed) < 0.25
    assert sum(d < 0.1 for d in fixed) / len(fixed) > 0.95


def test_methyl_groups_are_staggered_like_the_real_ones(ubiquitin):
    _, methyl = compare_with_experiment(ubiquitin)
    assert len(methyl) > 100
    assert median(methyl) < 0.3


def test_existing_hydrogens_are_kept(ubiquitin):
    residues = ubiquitin.residues()
    result = add_hydrogens(residues)
    assert result.added == 1                           # solo el HE2 de la His68 (trae HD1)
    for original, new in zip(residues, result.residues):
        for name, atom in original.atoms.items():
            assert new.atoms[name] == atom


def test_disulfide_cysteines_have_no_sh_hydrogen(hemoglobin):
    residues = hemoglobin.residues()
    cys = [i for i, r in enumerate(residues) if r.resn == "CYS"]
    a, b = residues[cys[0]], residues[cys[1]]
    sg = a.atoms["SG"]
    moved = replace(b.atoms["SG"], x=sg.x + 2.03, y=sg.y, z=sg.z)          # S-S a 2,03 Å
    residues[cys[1]] = replace(b, atoms={**b.atoms, "SG": moved})
    result = add_hydrogens(residues)
    assert "HG" not in result.residues[cys[0]].atoms
    assert "HG" not in result.residues[cys[1]].atoms


def test_chain_break_starts_a_new_n_terminus(hemoglobin):
    residues = hemoglobin.residues()
    del residues[50]                                   # falta un aminoácido: cadena cortada
    after = add_hydrogens(residues).residues[50]
    assert "H" not in after.atoms
    assert len({"H1", "H2", "H3"} & set(after.atoms)) == (2 if after.resn == "PRO" else 3)


def test_missing_reference_atoms_are_reported(hemoglobin):
    residues = hemoglobin.residues()
    broken = residues[5]
    residues[5] = replace(broken, atoms={n: a for n, a in broken.atoms.items() if n != "CB"})
    result = add_hydrogens(residues)
    assert any(resi == broken.resi for resi, _ in result.skipped)


# ---------- Molécula 3D (necesita NumPy) ----------

@pytest.fixture(scope="module")
def protein(hemoglobin):
    pytest.importorskip("numpy")
    from dnainvent.scene.protein import build_protein
    return build_protein(hemoglobin)


def test_protein_molecule_contents(protein):
    from collections import Counter
    assert protein.length == 147 and protein.sequence == HBB
    assert Counter(protein.elements) == {"H": 1127, "C": 729, "O": 202, "N": 196, "S": 4}
    assert protein.calculated_h.sum() == 1127
    assert protein.confidence is not None and protein.confidence.max() <= 100
    assert abs(protein.positions.mean(axis=0)).max() < 1e-3     # centrada
    assert set(protein.part) == {"backbone", "sidechain"}


def test_protein_bonds_are_chemically_sane(protein):
    import numpy as np
    lengths = np.linalg.norm(protein.positions[protein.bonds[:, 0]]
                             - protein.positions[protein.bonds[:, 1]], axis=1) * 10   # Å
    assert 0.9 < lengths.min() and lengths.max() < 1.9       # de O-H (0,96) a C-S (1,82)
    valence = np.bincount(protein.bonds.ravel(), minlength=protein.n_atoms)
    assert (valence > 0).all()                               # ningún átomo suelto
    hydrogens = np.array(protein.elements) == "H"
    assert (valence[hydrogens] == 1).all()                   # cada H, un enlace


def test_no_hydrogen_clashes(protein):
    import numpy as np
    p = protein.positions * 10
    bonded = {tuple(sorted(b)) for b in protein.bonds.tolist()}
    is_h = np.array(protein.elements) == "H"
    for i in np.where(is_h)[0]:
        d = np.linalg.norm(p - p[i], axis=1)
        close = [j for j in np.where(d < 1.6)[0] if j != i and tuple(sorted((i, j))) not in bonded]
        assert not close, (protein.names[i], int(protein.pair_index[i]))


def test_protein_window(hemoglobin):
    pytest.importorskip("numpy")
    from dnainvent.scene.protein import build_protein
    part = build_protein(hemoglobin, start=10, count=20)
    assert (part.start, part.end, part.length) == (10, 30, 147)
    assert part.sequence == HBB[10:30]
    assert set(part.pair_index.tolist()) == set(range(10, 30))
    assert part.residue_numbers == list(range(11, 31))       # numeración del archivo
    peptide = sum(1 for a, b in part.bonds.tolist() if {part.names[a], part.names[b]} == {"C", "N"}
                  and part.pair_index[a] != part.pair_index[b])
    assert peptide == 19


def test_experimental_structure_has_no_confidence(ubiquitin):
    pytest.importorskip("numpy")
    from dnainvent.scene.protein import build_protein
    molecule = build_protein(ubiquitin)
    assert molecule.confidence is None
    assert not molecule.calculated_h[:20].any()              # sus H vienen del archivo


def test_protein_molecule_has_the_molecule_interface(protein):
    pytest.importorskip("atomicworld")
    from dataclasses import fields
    from dnainvent.scene.molecule import Molecule
    assert {f.name for f in fields(Molecule)} <= {f.name for f in fields(type(protein))}


# ---------- Confianza ----------

@pytest.mark.parametrize("value, band", [(97.2, "very_high"), (90, "very_high"),
                                         (89.9, "confident"), (70, "confident"),
                                         (55, "low"), (49.9, "very_low"), (0, "very_low")])
def test_confidence_bands(value, band):
    assert confidence.band(value) == band


# ---------- Caché ----------

def test_cache_roundtrip_and_expiry(tmp_path):
    cache = Cache(tmp_path)
    assert cache.get("afdb", "X.json") is None
    cache.put("afdb", "X.json", b"hola")
    assert cache.get("afdb", "X.json") == b"hola"
    old = time.time() - 3600
    os.utime(cache.path("afdb", "X.json"), (old, old))
    assert cache.get("afdb", "X.json", max_age=60) is None
    assert cache.get("afdb", "X.json", max_age=7200) == b"hola"


def test_cache_keys_are_safe(tmp_path):
    cache = Cache(tmp_path)
    for namespace, key in (("../fuera", "a/b:c"), ("..", ".."), (".", "x"), ("a", "")):
        path = cache.path(namespace, key).resolve()
        assert tmp_path.resolve() in path.parents, (namespace, key)
        assert path.parent.parent == tmp_path.resolve()
    assert key_for("MKV") == key_for("MKV") != key_for("MKA")


def test_cache_folder_from_environment(monkeypatch, tmp_path):
    monkeypatch.setenv("DNAINVENT_CACHE", str(tmp_path))
    assert Cache().root == tmp_path


# ---------- AlphaFold DB (sin red) ----------

NEW_API = [  # campos de la v6, con un complejo (homodímero) que no se debe elegir
    {"modelEntityId": "AF-0000000065932696", "latestVersion": 1, "sequence": "MV",
     "pdbUrl": "https://example.org/dimer.pdb"},
    {"modelEntityId": "AF-P68871-F1", "latestVersion": 6, "sequence": HBB,
     "globalMetricValue": 97.19,
     "pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-P68871-F1-model_v6.pdb"},
]
OLD_API = [{"entryId": "AF-P68871-F1", "latestVersion": 4, "uniprotSequence": HBB}]


class FakeHttp:
    def __init__(self, routes):
        self.routes, self.calls = routes, []

    def __call__(self, url):
        self.calls.append(url)
        answer = self.routes.get(url)
        if isinstance(answer, Exception):
            raise answer
        if answer is None:
            raise urllib.error.HTTPError(url, 404, "Not Found", None, None)
        return answer


def api(accession="P68871"):
    return f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"


@pytest.fixture
def afdb_routes():
    return {api(): json.dumps(NEW_API).encode(),
            NEW_API[1]["pdbUrl"]: HEMOGLOBIN.read_bytes()}


def test_alphafold_entry_and_structure(tmp_path, afdb_routes):
    db = AlphaFoldDB(Cache(tmp_path), FakeHttp(afdb_routes))
    entry, structure = db.structure("p68871")
    assert (entry.model_id, entry.version, entry.mean_plddt) == ("AF-P68871-F1", 6, 97.19)
    assert entry.sequence == structure.sequence() == HBB


def test_alphafold_uses_the_cache(tmp_path, afdb_routes):
    http = FakeHttp(afdb_routes)
    AlphaFoldDB(Cache(tmp_path), http).structure("P68871")
    AlphaFoldDB(Cache(tmp_path), http).structure("P68871")
    assert len(http.calls) == 2                       # ficha + PDB, solo la primera vez


def test_alphafold_old_field_names_and_missing_pdb_url(tmp_path):
    url = "https://alphafold.ebi.ac.uk/files/AF-P68871-F1-model_v4.pdb"
    http = FakeHttp({api(): json.dumps(OLD_API).encode(), url: HEMOGLOBIN.read_bytes()})
    entry, structure = AlphaFoldDB(Cache(tmp_path), http).structure("P68871")
    assert entry.pdb_url == url and entry.version == 4 and entry.sequence == HBB
    assert len(structure.atoms) == 1131


@pytest.mark.parametrize("answer, code", [
    (None, "protein_not_found"),                                       # 404
    (b"[]", "protein_not_found"),                                      # lista vacía
    (b"<html>error</html>", "service_unavailable"),                    # no es JSON
    (urllib.error.HTTPError(api(), 500, "Error", None, None), "service_unavailable"),
    (urllib.error.URLError("sin red"), "service_unavailable"),
    (TimeoutError(), "service_unavailable"),
])
def test_alphafold_errors(tmp_path, answer, code):
    db = AlphaFoldDB(Cache(tmp_path), FakeHttp({api(): answer}))
    with pytest.raises(ServiceError) as info:
        db.entry("P68871")
    assert info.value.code == code


def test_invalid_accession_does_not_touch_the_network(tmp_path):
    http = FakeHttp({})
    with pytest.raises(InvalidParameterError) as info:
        AlphaFoldDB(Cache(tmp_path), http).entry("HBB")
    assert info.value.code == "invalid_accession"
    assert http.calls == []


@pytest.mark.parametrize("accession", ["P68871", "Q9Y6K9", "A0A023GPI8", " p68871 "])
def test_valid_accessions(accession):
    assert validate_accession(accession) == accession.strip().upper()
