"""Tests de las fichas de la opción 9 (dnainvent.database.cards), sin red.

Se construye un catálogo de HBB con el NCBI simulado de test_database y se le añaden a mano
los casos difíciles que salieron en el catálogo real: un registro de ClinVar con varias
enfermedades a la vez, fármacos sobre la causa y para los síntomas, y uno retirado.
"""

import pytest

from dnainvent.database import cards, get_or_create_source, open_database
from test_catalog import make_builder
from test_database import FakeNCBI


@pytest.fixture
def conn(tmp_path):
    conn = open_database(tmp_path / "catalog.db")
    make_builder(conn, FakeNCBI(), tmp_path).run(["HBB"])
    yield conn
    conn.close()


def ids(conn):
    gene = conn.execute("SELECT id FROM genes WHERE symbol = 'HBB'").fetchone()[0]
    hbs = conn.execute("SELECT id FROM variants WHERE hgvs_c = 'c.20A>T'").fetchone()[0]
    sickle = conn.execute("SELECT id FROM diseases WHERE name_en = 'Hb SS disease'").fetchone()[0]
    return gene, hbs, sickle


# ---------------------------------------------------------------- piezas

@pytest.mark.parametrize("significance, expected", [
    ("pathogenic", "pathogenic"), ("pathogenic; other", "pathogenic"),
    ("pathogenic/likely pathogenic", "pathogenic"), ("likely pathogenic, low penetrance", "pathogenic"),
    ("uncertain significance", "uncertain"),
    ("conflicting classifications of pathogenicity; other; risk factor", "uncertain"),
    ("likely benign", "benign"), ("benign/likely benign", "benign"),
    ("protective", "special"), ("risk factor", "special"), ("association", "special"),
    ("drug response", "special"),
    ("other", "unclassified"), ("not provided", "unclassified"), (None, "unclassified"),
])
def test_color(significance, expected):
    assert cards.color(significance) == expected


@pytest.mark.parametrize("status, expected", [
    ("practice guideline", (100, "guideline")), ("reviewed by expert panel", (75, "expert_panel")),
    ("criteria provided, multiple submitters, no conflicts", (50, "multiple")),
    ("criteria provided, single submitter", (25, "single")),
    ("criteria provided, conflicting classifications", (25, "conflicting")),
    ("no assertion criteria provided", (0, "no_criteria")), (None, (0, "no_criteria")),
])
def test_review_level(status, expected):
    assert cards.review_level(status) == expected


def test_frequency_and_inheritance():
    assert cards.frequency(0.01298) == {"value": 0.01298, "one_in": 77}       # HbS en gnomAD
    assert cards.frequency(None) is None and cards.frequency(0) is None
    assert cards.inheritance("Autosomal dominant|Autosomal recessive") == [
        "autosomal_dominant", "autosomal_recessive"]
    assert cards.inheritance("Something new") == ["Something new"]
    assert cards.inheritance(None) == []


def test_generic_disease_names():
    for name in ("HBB-related disorder", "G6PD-related disorders", "Hereditary cancer-predisposing syndrome",
                 "BRCA2-related cancer predisposition"):
        assert cards.is_generic(name), name
    assert not cards.is_generic("Hb SS disease")


def test_withdrawal_scope():
    def drug(chembl_id, withdrawn):
        return {"chembl_id": chembl_id, "withdrawn": withdrawn}
    assert cards.withdrawal(drug("CHEMBL4101807", 1)) == {"scope": "worldwide", "year": 2024}   # voxelotor
    assert cards.withdrawal(drug("CHEMBL1", 1)) == {"scope": "some"}         # ChEMBL no dice dónde
    assert cards.withdrawal(drug("CHEMBL1", 0)) is None


# ---------------------------------------------------------------- fichas

def test_gene_card(conn):
    gene, _, _ = ids(conn)
    card = cards.gene_card(conn, gene)
    assert (card["symbol"], card["chromosome"], card["reference_mrna"]) == ("HBB", "11", "NM_000518.5")
    assert sum(card["colors"].values()) == card["total"]
    assert card["colors"]["pathogenic"] >= 1 and card["colors"]["uncertain"] >= 1
    assert [v["common_name"] for v in card["known_variants"]] == ["HbS"]


def test_variant_card_hbs(conn):
    _, hbs, _ = ids(conn)
    card = cards.variant_card(conn, hbs)
    assert (card["common_name"], card["color"], card["dbsnp"]) == ("HbS", "pathogenic", "rs334")
    effect = card["change"]["effect"]
    assert (effect.ref_codon, effect.alt_codon, effect.kind) == ("GAG", "GTG", "missense")
    names = {d["name"]["en"]: d for d in card["diseases"]}
    assert names["Hb SS disease"]["effect"] == "causal" and names["Hb SS disease"]["review"] == 25


def test_variant_card_hides_records_with_several_diseases(conn):
    """El caso real RCV000477892: un registro, sin criterios, con varias enfermedades. No debe
    hacer que HbS «cause» la susceptibilidad a la malaria en la lista principal."""
    _, hbs, sickle = ids(conn)
    source = get_or_create_source(conn, "clinvar")
    malaria = conn.execute("SELECT id FROM diseases WHERE name_en = 'Malaria, susceptibility to'").fetchone()[0]
    thal = conn.execute("SELECT id FROM diseases WHERE name_en = 'beta Thalassemia'").fetchone()[0]
    with conn:
        conn.execute("DELETE FROM variant_diseases WHERE variant_id = ? AND disease_id = ?", (hbs, malaria))
        for disease in (malaria, thal):
            # aunque tuviera criterios, al ser de varias enfermedades no cuenta
            conn.execute("INSERT INTO variant_diseases (variant_id, disease_id, effect, clinical_significance, "
                         "review_status, rcv_accession, source_id) VALUES (?, ?, 'causal', 'pathogenic', "
                         "'criteria provided, single submitter', 'RCV000477892', ?)", (hbs, disease, source))
    card = cards.variant_card(conn, hbs)
    shown = {d["name"]["en"] for d in card["diseases"]}
    assert "Malaria, susceptibility to" not in shown and "Hb SS disease" in shown
    assert card["other_diseases"] >= 1


def test_low_review_records_are_not_in_the_main_list(conn):
    _, hbs, sickle = ids(conn)
    with conn:
        conn.execute("UPDATE variant_diseases SET review_status = 'no assertion criteria provided' "
                     "WHERE variant_id = ? AND disease_id = ?", (hbs, sickle))
    card = cards.variant_card(conn, hbs)
    assert "Hb SS disease" not in {d["name"]["en"] for d in card["diseases"]}
    gene, _, _ = ids(conn)
    assert "Hb SS disease" not in {d["name"]["en"] for d in cards.main_diseases(conn, gene)}


def add_drug(conn, name, chembl_id, disease_id, *, molecule_type="Small molecule", gene_id=None,
             action=None, target=None, withdrawn=0, phase=4, evidence="label", year=2000):
    source = get_or_create_source(conn, "chembl")
    with conn:
        drug = conn.execute("INSERT INTO drugs (name, chembl_id, max_phase, molecule_type, first_approval, "
                            "withdrawn, source_id) VALUES (?, ?, 4, ?, ?, ?, ?)",
                            (name, chembl_id, molecule_type, year, withdrawn, source)).lastrowid
        conn.execute("INSERT INTO drug_diseases (drug_id, disease_id, max_phase, evidence, source_id) "
                     "VALUES (?, ?, ?, ?, ?)", (drug, disease_id, phase, evidence, source))
        if target:
            conn.execute("INSERT INTO drug_mechanisms (drug_id, target_chembl_id, target_name, target_type, "
                         "action, mechanism) VALUES (?, 'CHEMBLT', ?, 'SINGLE PROTEIN', ?, ?)",
                         (drug, target, action, f"{target} {action}"))
        if gene_id:
            conn.execute("INSERT INTO drug_targets (drug_id, gene_id, action, mechanism, source_id) "
                         "VALUES (?, ?, ?, ?, ?)", (drug, gene_id, action, f"{target} {action}", source))
    return drug


def test_disease_card_groups_drugs(conn):
    gene, _, sickle = ids(conn)
    add_drug(conn, "VOXELOTOR", "CHEMBL4101807", sickle, gene_id=gene, target="Hemoglobin HbA",
             action="positive modulator", withdrawn=1, year=2019)
    add_drug(conn, "EXAGAMGLOGENE AUTOTEMCEL", "CHEMBL5095474", sickle, molecule_type="Gene",
             target="BCL11 transcription factor A", action="gene editing negative modulator", year=2023)
    add_drug(conn, "HYDROXYUREA", "CHEMBL467", sickle, target="Ribonucleoside-diphosphate reductase RR1",
             action="inhibitor", year=1967)
    add_drug(conn, "APIXABAN", "CHEMBL231779", sickle, phase=3, evidence="trials", year=2011)
    add_drug(conn, "OLDTRIAL", "CHEMBL9", sickle, phase=4, evidence="trials")     # fase 4 de ensayo: no aprobado

    card = cards.disease_card(conn, sickle)
    assert card["status"] == "approved"
    assert [d["name"] for d in card["approved_cause"]] == ["EXAGAMGLOGENE AUTOTEMCEL", "VOXELOTOR"]  # retirado al final
    assert [d["name"] for d in card["approved_symptoms"]] == ["HYDROXYUREA"]
    assert {d["name"] for d in card["trials"]} == {"APIXABAN", "OLDTRIAL"}
    assert card["approved_cause"][1]["withdrawal"] == {"scope": "worldwide", "year": 2024}
    assert card["genes"][0]["symbol"] == "HBB"


def test_drug_card(conn):
    gene, _, sickle = ids(conn)
    drug = add_drug(conn, "VOXELOTOR", "CHEMBL4101807", sickle, gene_id=gene, target="Hemoglobin HbA",
                    action="positive modulator", withdrawn=1, year=2019)
    with conn:      # la misma diana dos veces con la misma acción: se enseña una
        conn.execute("INSERT OR IGNORE INTO drug_mechanisms (drug_id, target_chembl_id, target_name, target_type, "
                     "action, mechanism) VALUES (?, 'CHEMBLT2', 'Hemoglobin HbA', 'PROTEIN COMPLEX', "
                     "'positive modulator', 'otra')", (drug,))
    card = cards.drug_card(conn, drug)
    assert card["structure"] == "5E83" and card["withdrawal"]["scope"] == "worldwide"
    assert [(m["target_name"], m["action"]) for m in card["mechanisms"]] == [("Hemoglobin HbA", "positive modulator")]
    assert card["genes"][0]["symbol"] == "HBB"
    assert card["diseases"][0]["name"]["en"] == "Hb SS disease" and card["diseases"][0]["approved"]


def test_cards_work_on_a_read_only_catalog(tmp_path, conn):
    """La web y la consola abren el catálogo solo para leer: la comparación no debe escribir."""
    import sqlite3
    path = tmp_path / "catalog.db"
    conn.commit()
    ro = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    ro.row_factory = sqlite3.Row
    _, hbs, _ = ids(ro)
    assert cards.variant_card(ro, hbs)["change"] is not None
    ro.close()


def test_featured_variants_start_with_named_and_well_reviewed(conn):
    gene, hbs, _ = ids(conn)
    featured = cards.featured_variants(conn, gene)
    assert featured[0]["id"] == hbs                                   # HbS, con nombre propio
    assert all(r["common_name"] or cards.color(r["clinical_significance"]) == "pathogenic" for r in featured)
    reviews = [cards.review_level(r["review_status"])[0] for r in featured if not r["common_name"]]
    assert reviews == sorted(reviews, reverse=True)


# ---------------------------------------------------------------- consola

@pytest.fixture
def console_on_catalog(tmp_path, conn, monkeypatch):
    """La consola con un catálogo de prueba, abierto solo para leer, y teclas simuladas."""
    from i18n import Translator
    from ui import console
    conn.commit()
    monkeypatch.setenv("DNAINVENT_CATALOG", str(tmp_path / "catalog.db"))
    monkeypatch.setattr(console, "_catalog", None)

    def run(keys, language="es"):
        monkeypatch.setattr(console, "tr", Translator(language))
        answers = iter(keys)
        monkeypatch.setattr("builtins.input", lambda prompt="": next(answers))
        console.db_catalog()
    return run


@pytest.mark.parametrize("language", ["es", "ca", "en"])
def test_console_goes_through_every_card(console_on_catalog, capsys, language):
    # gen 1 -> variantes (Intro: las destacadas) -> la 1 (HbS) -> detalles técnicos -> volver
    # -> enfermedades -> la 1 -> volver -> volver -> salir
    console_on_catalog(["1", "1", "", "1", "1", "0", "2", "1", "", "0", ""], language)
    out = capsys.readouterr().out
    assert "HbS" in out and "GAG" in out and "GTG" in out and "VCV000000101" in out
    assert "{" not in out, "queda algún hueco sin rellenar"          # p. ej. {label:<16}
    disclaimer = {"es": "no consejo médico", "ca": "no consell mèdic", "en": "not medical advice"}[language]
    assert disclaimer in out


def test_console_without_catalog(tmp_path, monkeypatch, capsys):
    from ui import console
    monkeypatch.setenv("DNAINVENT_CATALOG", str(tmp_path / "no_existe.db"))
    monkeypatch.setattr(console, "_catalog", None)
    console.db_catalog()
    assert "build_catalog" in capsys.readouterr().out


def test_console_card_of_a_nonsense_variant(conn, console_on_catalog, capsys):
    """HBB c.52A>T (p.Lys18Ter): AAG pasa a TAG, una parada. Cerraba la consola (KeyError '*')."""
    gene, hbs, _ = ids(conn)
    row = conn.execute("SELECT * FROM variants WHERE id = ?", (hbs,)).fetchone()
    with conn:
        conn.execute("INSERT INTO variants (gene_id, sequence_id, hgvs_c, hgvs_p, variant_type, "
                     "clinical_significance, review_status, source_id, clinvar_id) VALUES "
                     "(?, ?, 'c.52A>T', 'p.Lys18Ter', 'snv', 'pathogenic', 'reviewed by expert panel', ?, 999)",
                     (row["gene_id"], row["sequence_id"], row["source_id"]))
    console_on_catalog(["1", "1", "c.52A>T", "0", "0", ""])
    out = capsys.readouterr().out
    assert "AAG" in out and "TAG" in out and "Codón de parada (Ter)" in out
