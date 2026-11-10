"""Tests de fármacos (ChEMBL), UniProt, RCSB y ligandos dentro de proteínas, sin red.

Las respuestas de ChEMBL y UniProt imitan su formato (JSON con listas "targets", "mechanisms",
"molecules", "drug_indications", "drug_warnings" y "page_meta"), con IDs de prueba: no son los
reales. El caso es la anemia falciforme: el voxelotor actúa sobre la hemoglobina (sale desde el
gen) y la hidroxiurea sobre otra proteína (solo sale desde la enfermedad).
"""

import json
import math
import urllib.error
import urllib.parse

import pytest

from dnainvent import InvalidParameterError, ServiceError
from dnainvent.database import chembl, open_database, queries
from dnainvent.database.chembl import ChEMBLProvider, efo_to_curie, phase
from dnainvent.protein import ligands
from dnainvent.protein.cache import Cache
from dnainvent.protein.clients.rcsb import RCSB
from dnainvent.protein.clients.uniprot import UniProt
from dnainvent.protein.pdb import parse_pdb

HB = {"target_chembl_id": "CHEMBL_HB", "pref_name": "Hemoglobin", "target_type": "PROTEIN COMPLEX",
      "target_components": [{"accession": "P69905"}, {"accession": "P68871"}]}
RNR = {"target_chembl_id": "CHEMBL_RNR", "pref_name": "Ribonucleoside-diphosphate reductase",
       "target_type": "PROTEIN COMPLEX", "target_components": [{"accession": "P23921"}]}
MECHANISMS = [
    {"molecule_chembl_id": "CHEMBL_VOX", "parent_molecule_chembl_id": "CHEMBL_VOX",
     "target_chembl_id": "CHEMBL_HB", "action_type": "MODULATOR",
     "mechanism_of_action": "Hemoglobin S polymerisation inhibitor"},
    {"molecule_chembl_id": "CHEMBL_HU", "parent_molecule_chembl_id": "CHEMBL_HU",
     "target_chembl_id": "CHEMBL_RNR", "action_type": "INHIBITOR",
     "mechanism_of_action": "Ribonucleoside-diphosphate reductase inhibitor"},
]
MOLECULES = [
    {"molecule_chembl_id": "CHEMBL_VOX", "pref_name": "VOXELOTOR", "molecule_type": "Small molecule",
     "max_phase": "4.0", "first_approval": 2019,
     "molecule_structures": {"canonical_smiles": "Cc1ccc", "standard_inchi_key": "AAA-BBB"}},
    {"molecule_chembl_id": "CHEMBL_HU", "pref_name": "HYDROXYUREA", "molecule_type": "Small molecule",
     "max_phase": "4.0", "first_approval": 1967, "molecule_structures": None},
    {"molecule_chembl_id": "CHEMBL_TRIAL", "pref_name": "TRIALOMAB", "molecule_type": "Antibody",
     "max_phase": "2.0"},
    {"molecule_chembl_id": "CHEMBL_STATIN", "pref_name": "ATORVASTATIN", "molecule_type": "Small molecule",
     "max_phase": "4.0", "first_approval": 1996},
]
# Como en el ChEMBL real (comprobado en noviembre de 2026): el efo_id lleva dos puntos y cada
# indicación dice de dónde sale (indication_refs). La estatina es un fármaco aprobado para otra
# cosa y probado en un ensayo de fase 4 con pacientes de anemia falciforme: NO está aprobada
# para ella.
LABEL = [{"ref_type": "FDA", "ref_id": "label"}]
TRIAL = [{"ref_type": "ClinicalTrials", "ref_id": "NCT00000000"}]
INDICATIONS = [
    {"molecule_chembl_id": "CHEMBL_VOX", "parent_molecule_chembl_id": "CHEMBL_VOX",
     "efo_id": "MONDO:0011382", "efo_term": "sickle cell anemia", "mesh_id": "D000755",
     "mesh_heading": "Anemia, Sickle Cell", "max_phase_for_ind": "4.0", "indication_refs": LABEL},
    {"molecule_chembl_id": "CHEMBL_HU", "parent_molecule_chembl_id": "CHEMBL_HU",
     "efo_id": "MONDO:0011382", "efo_term": "sickle cell anemia", "mesh_id": "D000755",
     "mesh_heading": "Anemia, Sickle Cell", "max_phase_for_ind": "4.0", "indication_refs": LABEL + TRIAL},
    {"molecule_chembl_id": "CHEMBL_HU", "parent_molecule_chembl_id": "CHEMBL_HU",
     "efo_id": "EFO:0000339", "efo_term": "chronic myelogenous leukemia", "mesh_id": "D015464",
     "mesh_heading": "Leukemia, Myelogenous, Chronic", "max_phase_for_ind": "4.0", "indication_refs": LABEL},
    {"molecule_chembl_id": "CHEMBL_TRIAL", "parent_molecule_chembl_id": "CHEMBL_TRIAL",
     "efo_id": "EFO:TEST", "efo_term": "Beta thalassemia", "max_phase_for_ind": "2.0", "indication_refs": TRIAL},
    {"molecule_chembl_id": "CHEMBL_STATIN", "parent_molecule_chembl_id": "CHEMBL_STATIN",
     "efo_id": "MONDO:0011382", "efo_term": "sickle cell anemia", "mesh_id": "D000755",
     "mesh_heading": "Anemia, Sickle Cell", "max_phase_for_ind": "4.0", "indication_refs": TRIAL},
]
WARNINGS = [{"molecule_chembl_id": "CHEMBL_VOX", "warning_type": "Withdrawn", "warning_year": 2024}]


def page(key, items):
    return {key: items, "page_meta": {"total_count": len(items), "limit": 1000, "offset": 0}}


class FakeEBI:
    def __init__(self):
        self.calls = []

    def __call__(self, url: str) -> bytes:
        parsed = urllib.parse.urlparse(url)
        q = dict(urllib.parse.parse_qsl(parsed.query))
        self.calls.append(url)
        if parsed.netloc == "rest.uniprot.org":
            gene = q["query"].split()[0].split(":")[1]
            accession = {"HBB": "P68871", "HBA1": "P69905"}.get(gene.upper())   # UniProt no distingue mayúsculas
            results = [{"primaryAccession": accession, "genes": [{"geneName": {"value": gene}}]}] if accession else []
            return json.dumps({"results": results}).encode()
        resource = parsed.path.rsplit("/", 1)[-1].removesuffix(".json")

        def among(field):
            return set(q.get(f"{field}__in", "").split(","))

        if resource == "target":
            if "target_components__accession" in q:
                items = [t for t in (HB, RNR) if any(c["accession"] == q["target_components__accession"]
                                                     for c in t["target_components"])]
            else:
                items = [t for t in (HB, RNR) if t["target_chembl_id"] in among("target_chembl_id")]
            return json.dumps(page("targets", items)).encode()
        if resource == "mechanism":
            if "target_chembl_id__in" in q:
                items = [m for m in MECHANISMS if m["target_chembl_id"] in among("target_chembl_id")]
            else:
                items = [m for m in MECHANISMS if m["molecule_chembl_id"] in among("molecule_chembl_id")]
            return json.dumps(page("mechanisms", items)).encode()
        if resource == "molecule":
            return json.dumps(page("molecules", [m for m in MOLECULES
                                                 if m["molecule_chembl_id"] in among("molecule_chembl_id")])).encode()
        if resource == "drug_indication":
            if "efo_id" in q:
                items = [i for i in INDICATIONS if i["efo_id"] == q["efo_id"]]
            elif "mesh_id" in q:
                items = [i for i in INDICATIONS if i.get("mesh_id") == q["mesh_id"]]
            else:
                items = [i for i in INDICATIONS if i["molecule_chembl_id"] in among("molecule_chembl_id")]
            return json.dumps(page("drug_indications", items)).encode()
        if resource == "drug_warning":
            return json.dumps(page("drug_warnings", [w for w in WARNINGS
                                                     if w["molecule_chembl_id"] in among("molecule_chembl_id")])).encode()
        raise AssertionError(f"Petición inesperada: {url}")


@pytest.fixture
def conn():
    connection = open_database(":memory:")
    source = connection.execute("INSERT INTO sources (code, name, license) VALUES ('t', 't', 't')").lastrowid
    connection.execute("INSERT INTO organisms (taxon_id, scientific_name) VALUES (9606, 'Homo sapiens')")
    for symbol in ("HBB", "HBA1"):
        connection.execute("INSERT INTO genes (taxon_id, symbol, source_id) VALUES (9606, ?, ?)", (symbol, source))
    connection.execute("INSERT INTO diseases (mondo_id, name_en) VALUES ('MONDO:0011382', 'Hb SS disease')")
    connection.execute("INSERT INTO diseases (name_en) VALUES ('beta Thalassemia')")
    connection.commit()
    yield connection
    connection.close()


@pytest.fixture
def fake():
    return FakeEBI()


@pytest.fixture
def provider(fake, tmp_path):
    return ChEMBLProvider(http=fake, cache=Cache(tmp_path), sleep=lambda s: None)


def gene_id(conn, symbol):
    return conn.execute("SELECT id FROM genes WHERE symbol = ?", (symbol,)).fetchone()[0]


SICKLE = 1


# ---------------------------------------------------------------- utilidades

@pytest.mark.parametrize("value, expected", [("4.0", 4), (4, 4), ("0.5", 0.5), ("-1.0", None),
                                             (None, None), ("x", None)])
def test_phase(value, expected):
    assert phase(value) == expected


def test_efo_ids():
    assert efo_to_curie("MONDO_0011382") == "MONDO:0011382"
    assert efo_to_curie("MONDO:0011382") == "MONDO:0011382"
    assert chembl.curie_to_efo("MONDO:0011382") == "MONDO_0011382"


@pytest.mark.parametrize("refs, expected", [
    ([{"ref_type": "FDA"}], "label"), ([{"ref_type": "ClinicalTrials"}, {"ref_type": "DailyMed"}], "label"),
    ([{"ref_type": "ClinicalTrials"}], "trials"), ([{"ref_type": "PubMed"}], "other"), ([], None), (None, None),
    # USAN es el registro del NOMBRE del fármaco en EE. UU., no una aprobación
    ([{"ref_type": "USAN"}, {"ref_type": "ClinicalTrials"}], "trials"),
])
def test_evidence(refs, expected):
    assert chembl.evidence({"indication_refs": refs}) == expected


def test_uniprot_prefers_primary_gene_name(fake, tmp_path):
    assert UniProt(cache=Cache(tmp_path), http=fake).accession("hbb") == "P68871"
    assert UniProt(cache=Cache(tmp_path), http=fake).accession("NOEXISTE") is None


# ---------------------------------------------------------------- ChEMBL

def test_drugs_from_gene(conn, provider):
    stats = provider.import_gene_drugs(conn, gene_id(conn, "HBB"))
    assert (stats.accession, stats.targets, stats.drugs, stats.approved, stats.withdrawn) == ("P68871", 1, 1, 1, 1)
    assert conn.execute("SELECT uniprot_id FROM genes WHERE symbol = 'HBB'").fetchone()[0] == "P68871"
    [drug] = queries.gene_drugs(conn, gene_id(conn, "HBB"))
    assert (drug["name"], drug["action"], drug["withdrawn"], drug["first_approval"]) == ("VOXELOTOR", "modulator", 1, 2019)
    assert drug["mechanism"] == "Hemoglobin S polymerisation inhibitor"
    # Su indicación, enlazada a la enfermedad guardada por el ID de MONDO
    assert [d["name"] for d in queries.disease_drugs(conn, SICKLE)] == ["VOXELOTOR"]
    assert queries.disease_status(conn, SICKLE) == "withdrawn"


def test_complex_links_every_saved_gene(conn, provider):
    provider.import_gene_drugs(conn, gene_id(conn, "HBA1"))
    provider.import_gene_drugs(conn, gene_id(conn, "HBB"))
    # El complejo hemoglobina contiene las dos cadenas: el voxelotor queda en los dos genes
    assert [d["name"] for d in queries.gene_drugs(conn, gene_id(conn, "HBA1"))] == ["VOXELOTOR"]
    assert [d["name"] for d in queries.gene_drugs(conn, gene_id(conn, "HBB"))] == ["VOXELOTOR"]


def test_drugs_from_disease_find_other_targets(conn, provider):
    stats = provider.import_disease_drugs(conn, SICKLE)
    assert (stats.drugs, stats.approved) == (3, 3)
    drugs = {d["name"]: d for d in queries.disease_drugs(conn, SICKLE)}
    assert set(drugs) == {"HYDROXYUREA", "VOXELOTOR", "ATORVASTATIN"}
    assert (drugs["HYDROXYUREA"]["evidence"], drugs["ATORVASTATIN"]["evidence"]) == ("label", "trials")
    # Primero los aprobados para ella (y no retirados); la estatina, aunque sea fase 4, después
    assert [d["name"] for d in queries.disease_drugs(conn, SICKLE)] == ["HYDROXYUREA", "VOXELOTOR", "ATORVASTATIN"]
    assert drugs["HYDROXYUREA"]["mechanisms"] == "Ribonucleoside-diphosphate reductase inhibitor"
    assert drugs["HYDROXYUREA"]["phase_for_disease"] == 4
    assert queries.disease_status(conn, SICKLE) == "approved"      # la hidroxiurea no se retiró
    # La leucemia no está guardada: no se crea, se cuenta
    assert stats.indications_other == 1
    assert conn.execute("SELECT COUNT(*) FROM diseases").fetchone()[0] == 2
    # La hidroxiurea no actúa sobre HBB: no aparece en sus fármacos
    assert queries.gene_drugs(conn, gene_id(conn, "HBB")) == []


def test_phase_4_trials_alone_are_not_an_approval(conn, provider):
    """Solo la estatina (fase 4 en un ensayo, aprobada para otra cosa): no hay tratamiento aprobado."""
    with conn:
        source = conn.execute("SELECT id FROM sources LIMIT 1").fetchone()[0]
        statin = chembl._upsert_drug(conn, MOLECULES[3], False, source)
        conn.execute("INSERT INTO drug_diseases (drug_id, disease_id, max_phase, source_id, evidence) "
                     "VALUES (?, ?, 4, ?, 'trials')", (statin, SICKLE, source))
    assert queries.disease_status(conn, SICKLE) == "trials"


def test_mondo_with_underscore_is_tried_too(conn, provider, monkeypatch):
    asked = []
    monkeypatch.setattr(provider, "all", lambda r, k, params: asked.append(params.get("efo_id")) or [])
    provider.indications(efo_id="MONDO:0011382")
    assert asked == ["MONDO:0011382", "MONDO_0011382"]


def test_one_failing_way_of_searching_does_not_stop_the_rest(provider, monkeypatch):
    def flaky(resource, key, params):
        if params.get("efo_id") == "MONDO:0011382":
            raise ServiceError("service_unavailable", service="ChEMBL")
        return [{"molecule_chembl_id": "X", "via": params}]

    monkeypatch.setattr(provider, "all", flaky)
    found = provider.indications(efo_id="MONDO:0011382", mesh_id="D000755")
    assert [f["via"] for f in found] == [{"efo_id": "MONDO_0011382"}, {"mesh_id": "D000755"}]

    monkeypatch.setattr(provider, "all", lambda r, k, p: (_ for _ in ()).throw(
        ServiceError("service_unavailable", service="ChEMBL")))
    with pytest.raises(ServiceError):
        provider.indications(efo_id="MONDO:0011382")


def test_disease_status_trials_and_none(conn, provider):
    assert queries.disease_status(conn, 2) == "none"
    provider.molecules(["CHEMBL_TRIAL"])
    with conn:
        source = conn.execute("SELECT id FROM sources LIMIT 1").fetchone()[0]
        drug = chembl._upsert_drug(conn, MOLECULES[2], False, source)
        conn.execute("INSERT INTO drug_diseases (drug_id, disease_id, max_phase, source_id) VALUES (?, 2, 2, ?)",
                     (drug, source))
    assert queries.disease_status(conn, 2) == "trials"


def test_reimport_does_not_duplicate_and_uses_cache(conn, provider, fake):
    provider.import_disease_drugs(conn, SICKLE)
    calls = len(fake.calls)
    provider.import_disease_drugs(conn, SICKLE)
    assert len(fake.calls) == calls                       # todo de la caché
    assert conn.execute("SELECT COUNT(*) FROM drugs").fetchone()[0] == 3
    assert conn.execute("SELECT COUNT(*) FROM drug_mechanisms").fetchone()[0] == 2


def test_search_finds_drugs(conn, provider):
    provider.import_disease_drugs(conn, SICKLE)
    assert [(r["kind"], r["label"]) for r in queries.search(conn, "voxel")] == [("drug", "VOXELOTOR")]


def test_paging(provider, monkeypatch):
    pages = {0: {"targets": [{"target_chembl_id": "A"}, {"target_chembl_id": "B"}], "page_meta": {"total_count": 3}},
             2: {"targets": [{"target_chembl_id": "C"}], "page_meta": {"total_count": 3}}}
    monkeypatch.setattr(provider, "get", lambda resource, params: pages[params["offset"]])
    assert [t["target_chembl_id"] for t in provider.all("target", "targets", {})] == ["A", "B", "C"]


def test_in_filter_in_batches(provider, monkeypatch):
    seen = []
    monkeypatch.setattr(chembl, "IN_BATCH", 2)
    monkeypatch.setattr(provider, "all", lambda r, k, params: seen.append(params) or [])
    provider.molecules(["C", "A", "B", None])
    assert [p["molecule_chembl_id__in"] for p in seen] == ["A,B", "C"]


def test_slow_chembl_is_retried_and_announced(conn, tmp_path, fake):
    """Como el 10-11-2026: errores 500 y tiempos agotados antes de responder."""
    failures = [urllib.error.HTTPError("u", 500, "Error", None, None), TimeoutError("lento")]
    waits, notices = [], []

    def flaky(url):
        if failures:
            raise failures.pop(0)
        return fake(url)

    provider = ChEMBLProvider(http=flaky, cache=Cache(tmp_path), sleep=waits.append)
    provider.import_disease_drugs(conn, SICKLE, progress=lambda key, data: notices.append((key, data)))
    retries = [data for key, data in notices if key == "database.stage.chembl_retry"]
    assert [r["attempt"] for r in retries] == [2, 3]
    assert [w for w in waits if w >= 5] == [5, 10]
    assert {d["name"] for d in queries.disease_drugs(conn, SICKLE)} >= {"HYDROXYUREA"}


def test_chembl_down(conn, tmp_path):
    def down(url):
        raise urllib.error.URLError("sin red")

    provider = ChEMBLProvider(http=down, cache=Cache(tmp_path), sleep=lambda s: None)
    with pytest.raises(ServiceError) as info:
        provider.import_disease_drugs(conn, SICKLE)
    assert info.value.code == "service_unavailable"


# ---------------------------------------------------------------- PDB y ligandos

def atom_line(record, serial, name, resn, chain, resi, x, y, z, element):
    name = f" {name:<3}" if len(name) < 4 else name
    return (f"{record:<6}{serial:5d} {name} {resn:>3} {chain}{resi:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}")


def complex_pdb() -> str:
    """Una «proteína» de 2 aminoácidos, un fármaco de 6 átomos (LIG) unido al N del primero,
    un hemo, un sulfato y una molécula de agua."""
    lines = ["TITLE     PRUEBA CON FARMACO",
             "DBREF  1ABC A    1   141  UNP    P69905   HBA_HUMAN        1    141"]
    n = 0

    def add(*args):
        nonlocal n
        n += 1
        lines.append(atom_line(args[0], n, *args[1:]))

    add("ATOM", "N", "VAL", "A", 1, 0.0, 0.0, 0.0, "N")
    add("ATOM", "CA", "VAL", "A", 1, 1.5, 0.0, 0.0, "C")
    add("ATOM", "N", "LEU", "A", 2, 0.0, 6.0, 0.0, "N")
    add("ATOM", "CD1", "LEU", "A", 2, 4.0, 6.0, 0.0, "C")
    # El fármaco: C1 a 1,4 Å del N de la Val1 (enlace), O2 a 3,0 Å del N de la Leu2 (puente H),
    # C3 a 3,8 Å del CD1 de la Leu2 (hidrófobo)
    add("HETATM", "C1", "LIG", "A", 201, -1.4, 0.0, 0.0, "C")
    add("HETATM", "O2", "LIG", "A", 201, 0.0, 3.0, 0.0, "O")
    add("HETATM", "C3", "LIG", "A", 201, 4.0, 2.2, 0.0, "C")
    for i in range(3):
        add("HETATM", f"C{i + 4}", "LIG", "A", 201, -3.0 - i, -3.0, 0.0, "C")
    add("HETATM", "H1", "LIG", "A", 201, -1.4, 1.0, 0.0, "H")
    for i in range(10):
        add("HETATM", f"C{i}", "HEM", "A", 142, 20.0 + i, 20.0, 0.0, "C")
    add("HETATM", "S", "SO4", "A", 301, 30.0, 0.0, 0.0, "S")
    add("HETATM", "O", "HOH", "A", 401, 0.5, 0.5, 0.5, "O")
    lines.append("END")
    return "\n".join(lines) + "\n"


def test_dbref_chain_accessions():
    assert parse_pdb(complex_pdb()).chain_accessions == {"A": "P69905"}


def test_find_ligands():
    structure = parse_pdb(complex_pdb())
    assert [(lig.code, lig.kind) for lig in ligands.find_ligands(structure)] == [("LIG", "drug"), ("HEM", "cofactor")]
    everything = ligands.find_ligands(structure, include_additives=True)
    assert [lig.code for lig in everything] == ["LIG", "HEM", "SO4"]      # el agua nunca
    drug = everything[0]
    assert (drug.key, len(drug.heavy_atoms)) == ("LIG/A201", 6)


def test_contacts_and_pocket():
    structure = parse_pdb(complex_pdb())
    drug = ligands.find_ligands(structure)[0]
    found = ligands.contacts(structure, drug)
    assert all(c.protein_atom.element != "H" and not c.protein_atom.hetero for c in found)
    first = found[0]
    assert (first.kind, first.protein_atom.name, first.distance) == ("covalent", "N", 1.4)
    pocket = ligands.pocket(structure, drug)
    assert [(r.resn, r.resi, r.kind) for r in pocket] == [("VAL", 1, "covalent"), ("LEU", 2, "hbond")]
    leu = pocket[1]
    assert leu.distance == 3.0 and leu.contacts >= 2         # puente H con el N y hidrófobo con CD1
    kinds = {(c.ligand_atom.name, c.protein_atom.name): c.kind for c in found}
    assert kinds[("C3", "CD1")] == "hydrophobic"
    assert math.isclose(min(c.distance for c in found if c.kind == "hbond"), 3.0)


def test_rcsb_downloads_once_and_validates(tmp_path):
    calls = []

    def http(url):
        calls.append(url)
        return complex_pdb().encode()

    pdb = RCSB(cache=Cache(tmp_path), http=http)
    assert pdb.structure("5e83").chain_accessions == {"A": "P69905"}
    pdb.structure("5E83")
    assert calls == ["https://files.rcsb.org/download/5E83.pdb"]
    with pytest.raises(InvalidParameterError):
        pdb.structure("../etc")

    def missing(url):
        raise urllib.error.HTTPError(url, 404, "Not Found", None, None)

    with pytest.raises(ServiceError) as info:
        RCSB(cache=Cache(tmp_path), http=missing).structure("9ZZZ")
    assert info.value.code == "structure_not_found"
