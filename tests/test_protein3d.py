"""Tests de la mutación en la proteína 3D (dnainvent.database.protein3d) y de los enlaces de las
fichas, sin red: la estructura de AlphaFold es la de tests/data (hemoglobina beta)."""

from pathlib import Path

import pytest

from dnainvent.database import cards, protein3d
from dnainvent.errors import ServiceError
from dnainvent.protein.clients.afdb import AlphaFoldEntry
from dnainvent.protein.pdb import read_pdb

DATA = Path(__file__).parent / "data"


def hbb():
    return read_pdb(DATA / "AF-P68871-F1-model_v6.pdb")


class FakeAlphaFold:
    def __init__(self, fail=False):
        self.fail, self.asked = fail, []

    def structure(self, accession):
        self.asked.append(accession)
        if self.fail:
            raise ServiceError("service_unavailable", service="AlphaFold DB")
        s = hbb()
        return AlphaFoldEntry(accession, f"AF-{accession}-F1", 6, s.sequence(), "", 95.0), s


@pytest.mark.parametrize("hgvs_p, expected", [
    ("p.Glu7Val", (7, "E", "V", "missense")),
    ("p.Lys18Ter", (18, "K", "*", "nonsense")),
    ("p.Phe508del", (508, "F", None, "deletion")),
    ("p.Glu23fs", (23, "E", None, "frameshift")),
    ("p.(Ser1982ArgfsTer22)", (1982, "S", None, "frameshift")),
    ("p.Ter148Glu", (148, "*", "E", "stop_lost")),
])
def test_parse_hgvs_p(hgvs_p, expected):
    change = protein3d.parse_hgvs_p(hgvs_p)
    assert (change["position"], change["ref"], change["alt"], change["kind"]) == expected
    assert "(" not in change["label"]


@pytest.mark.parametrize("hgvs_p", ["p.?", "p.Leu4=", None, "", "c.20A>T"])
def test_parse_hgvs_p_without_an_amino_acid(hgvs_p):
    assert protein3d.parse_hgvs_p(hgvs_p) is None


def test_change_properties():
    # HbS: el glutamato (negativo, se lleva bien con el agua) pasa a valina (la repele)
    assert protein3d.change_properties("E", "V") == ["loses_negative", "polar_to_hydrophobic"]
    assert protein3d.change_properties("E", "K") == ["charge_reversed"]          # HbC
    assert protein3d.change_properties("G", "W") == ["polar_to_hydrophobic", "bigger", "from_glycine"]
    assert "to_proline" in protein3d.change_properties("L", "P")
    assert protein3d.change_properties("R", "C")[-1] == "cysteine"
    assert protein3d.change_properties("K", "*") == []


def test_analyze_hbs_is_on_the_surface():
    site = protein3d.analyze(hbb(), hbb().sequence(),
                             {"position": 7, "ref": "E", "alt": "V", "kind": "missense", "label": "Glu7Val"})
    assert site["available"] and (site["index"], site["number"], site["length"]) == (6, 7, 147)
    assert site["burial"] == "exposed" and site["confidence"] == "very_high"
    assert site["changes"] == ["loses_negative", "polar_to_hydrophobic"]
    assert site["contacts"] and all(abs(c["index"] - 6) > 1 for c in site["contacts"])
    assert (site["start"], site["count"]) == (0, 147)         # cabe entera en el microscopio


def test_analyze_nonsense_counts_what_is_lost():
    site = protein3d.analyze(hbb(), hbb().sequence(),
                             {"position": 18, "ref": "K", "alt": "*", "kind": "nonsense", "label": "Lys18Ter"})
    assert site["lost"] == 130


def test_analyze_refuses_the_wrong_amino_acid():
    """Si la estructura no tiene el aminoácido que dice la variante, no se marca nada."""
    site = protein3d.analyze(hbb(), hbb().sequence(),
                             {"position": 7, "ref": "W", "alt": "V", "kind": "missense", "label": "Trp7Val"})
    assert site == {"available": False, "reason": "mismatch"}
    far = protein3d.analyze(hbb(), hbb().sequence(),
                            {"position": 400, "ref": "E", "alt": "V", "kind": "missense", "label": "x"})
    assert far["reason"] == "not_in_structure"


def test_map_position_aligns_different_versions():
    refseq = hbb().sequence()
    shorter = refseq[1:]                       # la estructura sin la metionina inicial
    assert protein3d.map_position(refseq, refseq, 7) == 6
    assert protein3d.map_position(refseq, shorter, 7) == 5
    assert protein3d.map_position(refseq, refseq, 0) is None


def test_long_window_is_centered_on_the_site():
    structure = hbb()
    site = protein3d.analyze(structure, structure.sequence(),
                             {"position": 140, "ref": structure.sequence()[139], "alt": "A",
                              "kind": "missense", "label": "x"})
    assert site["start"] + site["count"] <= site["length"] and site["start"] <= 139


# ---------------------------------------------------------------- con el catálogo
@pytest.fixture
def conn(tmp_path):
    from test_catalog import make_builder
    from test_database import FakeNCBI
    from dnainvent.database import open_database
    conn = open_database(tmp_path / "catalog.db")
    make_builder(conn, FakeNCBI(), tmp_path).run(["HBB"])
    yield conn
    conn.close()


def with_protein(conn):
    """El gen con su accesión de UniProt y su proteína de RefSeq (la de la estructura)."""
    with conn:
        conn.execute("UPDATE genes SET uniprot_id = 'P68871' WHERE symbol = 'HBB'")
        gene = conn.execute("SELECT id FROM genes WHERE symbol = 'HBB'").fetchone()[0]
        updated = conn.execute("UPDATE sequences SET residues = ? WHERE gene_id = ? AND mol_type = 'protein'",
                               (hbb().sequence(), gene)).rowcount
        if not updated:
            source = conn.execute("SELECT id FROM sources LIMIT 1").fetchone()[0]
            conn.execute("INSERT INTO sequences (gene_id, accession, version, mol_type, residues, source_id) "
                         "VALUES (?, 'NP_000509', 1, 'protein', ?, ?)", (gene, hbb().sequence(), source))


def hbs_id(conn):
    return conn.execute("SELECT id FROM variants WHERE hgvs_c = 'c.20A>T'").fetchone()[0]


def test_variant_site_on_the_catalog(conn):
    afdb = FakeAlphaFold()
    with_protein(conn)
    site = protein3d.variant_site(conn, hbs_id(conn), afdb)
    assert site["available"] and site["label"] == "Glu7Val" and afdb.asked == ["P68871"]
    assert site["model_id"] == "AF-P68871-F1" and site["structure"] is not None

    failing = protein3d.variant_site(conn, hbs_id(conn), FakeAlphaFold(fail=True))
    assert (failing["available"], failing["reason"]) == (False, "service")

    with conn:            # una proteína de más de 2.700 aminoácidos: ni se pide a AlphaFold
        conn.execute("UPDATE sequences SET residues = ? WHERE mol_type = 'protein'", ("M" * 3000,))
    afdb = FakeAlphaFold()
    site = protein3d.variant_site(conn, hbs_id(conn), afdb)
    assert (site["reason"], site["length"], afdb.asked) == ("too_long", 3000, [])


def test_console_shows_the_site(conn, tmp_path, monkeypatch, capsys):
    from i18n import Translator
    from ui import console
    with_protein(conn)
    monkeypatch.setenv("DNAINVENT_CATALOG", str(tmp_path / "catalog.db"))
    monkeypatch.setattr(console, "_catalog", None)
    monkeypatch.setattr(console, "AlphaFoldDB", FakeAlphaFold)
    for language, words in (("es", "en la superficie"), ("ca", "a la superfície"), ("en", "on the surface")):
        monkeypatch.setattr(console, "tr", Translator(language))
        answers = iter(["3", "0", "0"])
        monkeypatch.setattr("builtins.input", lambda prompt="": next(answers))
        console.show_variant_card(hbs_id(conn))
        out = capsys.readouterr().out
        assert "Glu7Val" in out and words in out and "AF-P68871-F1" in out
        assert "{" not in out.replace("{}", "")                 # ningún texto sin rellenar


def test_links():
    gene = cards.gene_links({"symbol": "HBB", "ncbi_gene_id": 3043, "uniprot_id": "P68871"})
    assert ("NCBI Gene", "https://www.ncbi.nlm.nih.gov/gene/3043") in gene
    assert ("AlphaFold DB", "https://alphafold.ebi.ac.uk/entry/P68871") in gene
    variant = cards.variant_links({"clinvar_id": 15333, "dbsnp": "rs334"})
    assert variant == [("ClinVar", "https://www.ncbi.nlm.nih.gov/clinvar/variation/15333/"),
                       ("dbSNP", "https://www.ncbi.nlm.nih.gov/snp/rs334")]
    disease = cards.disease_links({"mondo_id": "MONDO:0011382", "orpha_codes": ["232"]})
    assert ("Orphanet", "https://www.orpha.net/en/disease/detail/232") in disease
    drug = cards.drug_links({"chembl_id": "CHEMBL4101807", "structure": "5E83"})
    assert drug[1] == ("RCSB PDB", "https://www.rcsb.org/structure/5E83")
    assert cards.variant_links({}) == []
