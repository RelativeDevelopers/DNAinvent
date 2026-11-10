"""Tests de la base de datos de genes reales (dnainvent/database), sin red.

Las respuestas falsas imitan el formato de las E-utilities del NCBI. El ARNm de HBB está
ACORTADO (no es el registro real), pero los primeros codones sí son los reales, así que la
variante de la anemia falciforme cae donde debe. tests/data/clinvar_esummary_hbb.json sí es
una respuesta real de ClinVar (3 variantes de HBB, noviembre de 2026).
"""

import json
import sqlite3
import urllib.error
import urllib.parse
from pathlib import Path

import pytest

from dnainvent import InvalidParameterError, ServiceError
from dnainvent.database import SCHEMA_PATH, clinvar, open_database, queries
from dnainvent.database.clinvar import (ClinVarProvider, classify_effect, parse_summary,
                                        parse_variation_name, parse_vcv)
from dnainvent.database.ncbi import NCBIProvider, parse_insdseq, parse_location
from dnainvent.database.variants import locate, side_by_side

DATA = Path(__file__).parent / "data"
REAL_SAMPLE = json.loads((DATA / "clinvar_esummary_hbb.json").read_text(encoding="utf-8"))["result"]

UTR5 = "ACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACC"                      # 50 nt
CDS = ("ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAG"
       "GTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTAA")
UTR3 = "GCTCGCTTTCTTGCTGTCCAATTTCTATTAAAGG"
MRNA = UTR5 + CDS + UTR3
CDS_START, CDS_END = len(UTR5) + 1, len(UTR5) + len(CDS)
PROTEIN = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGR"


def insdseq(accession="NM_000518.5", mrna=MRNA, prefix="INSD") -> bytes:
    p = prefix
    return f"""<?xml version="1.0"?>
<{p}Set><{p}Seq>
  <{p}Seq_moltype>mRNA</{p}Seq_moltype>
  <{p}Seq_definition>Homo sapiens hemoglobin subunit beta (HBB), mRNA</{p}Seq_definition>
  <{p}Seq_accession-version>{accession}</{p}Seq_accession-version>
  <{p}Seq_feature-table>
    <{p}Feature><{p}Feature_key>source</{p}Feature_key>
      <{p}Feature_location>1..{len(mrna)}</{p}Feature_location>
      <{p}Feature_quals><{p}Qualifier><{p}Qualifier_name>db_xref</{p}Qualifier_name>
        <{p}Qualifier_value>taxon:9606</{p}Qualifier_value></{p}Qualifier></{p}Feature_quals>
    </{p}Feature>
    <{p}Feature><{p}Feature_key>CDS</{p}Feature_key>
      <{p}Feature_location>{CDS_START}..{CDS_END}</{p}Feature_location>
      <{p}Feature_quals>
        <{p}Qualifier><{p}Qualifier_name>product</{p}Qualifier_name>
          <{p}Qualifier_value>hemoglobin subunit beta</{p}Qualifier_value></{p}Qualifier>
        <{p}Qualifier><{p}Qualifier_name>protein_id</{p}Qualifier_name>
          <{p}Qualifier_value>NP_000509.1</{p}Qualifier_value></{p}Qualifier>
        <{p}Qualifier><{p}Qualifier_name>translation</{p}Qualifier_name>
          <{p}Qualifier_value>{PROTEIN}</{p}Qualifier_value></{p}Qualifier>
      </{p}Feature_quals>
    </{p}Feature>
  </{p}Seq_feature-table>
  <{p}Seq_sequence>{mrna.lower()}</{p}Seq_sequence>
</{p}Seq></{p}Set>""".encode()


GENE_SUMMARY = {"result": {"uids": ["3043"], "3043": {
    "name": "HBB", "description": "hemoglobin subunit beta", "chromosome": "11",
    "organism": {"scientificname": "Homo sapiens", "commonname": "human", "taxid": 9606},
    "genomicinfo": [{"chraccver": "NC_000011.10", "chrstart": 5227070, "chrstop": 5225463}]}}}


def cv_doc(name, obj_type, significance, traits=(), rs=None, somatic=None):
    """Un documento esummary de ClinVar (los IDs de MedGen y MONDO son de prueba)."""
    return {
        "title": name, "obj_type": obj_type,
        "variation_set": [{"variation_name": name,
                           "variation_xrefs": [{"db_source": "dbSNP", "db_id": rs}] if rs else []}],
        "germline_classification": {"description": significance, "review_status": "reviewed",
                                    "trait_set": [{"trait_name": t, "trait_xrefs": xrefs}
                                                  for t, xrefs in traits]},
        "oncogenicity_classification": {"description": somatic or ""},
    }


SICKLE = ("Hb SS disease", [{"db_source": "MedGen", "db_id": "C1"},
                            {"db_source": "MONDO", "db_id": "MONDO:T1"}])
MALARIA = ("Malaria, susceptibility to", [{"db_source": "MedGen", "db_id": "C2"}])
THAL = ("beta Thalassemia", [{"db_source": "MedGen", "db_id": "C3"},
                             {"db_source": "MONDO", "db_id": "MONDO:T3"}])

CLINVAR_DOCS = {
    # HbS: patogénica en global, pero protectora frente a la malaria
    "101": cv_doc("NM_000518.5(HBB):c.20A>T (p.Glu7Val)", "single nucleotide variant",
                  "Pathogenic; other", [SICKLE, MALARIA], rs="334"),
    "102": cv_doc("NM_000518.5(HBB):c.20A>G (p.Glu7Gly)", "single nucleotide variant",
                  "Likely benign", [("not provided", [])]),
    "103": cv_doc("NM_000518.4(HBB):c.79G>A (p.Glu27Lys)", "single nucleotide variant",
                  "Pathogenic/Likely pathogenic", [THAL, SICKLE]),
    "104": cv_doc("GRCh38/hg38 11p15.4(chr11:5200000-5300000)x1", "copy number loss",
                  "Pathogenic", [THAL]),
    "105": cv_doc("NM_000518.5(HBB):c.93-21G>A", "single nucleotide variant",
                  "Conflicting classifications of pathogenicity", [THAL]),
    "106": cv_doc("NM_000519.4(HBD):c.-127T>C", "single nucleotide variant",
                  "Benign", [("not provided", [])]),
}


def rcv(accession, condition, medgen, significance):
    return f"""<RCVAccession Accession="{accession}" Version="1">
      <ClassifiedConditionList TraitSetID="1">
        <ClassifiedCondition DB="MedGen" ID="{medgen}">{condition}</ClassifiedCondition>
      </ClassifiedConditionList>
      <RCVClassifications><GermlineClassification>
        <ReviewStatus>criteria provided, single submitter</ReviewStatus>
        <Description DateLastEvaluated="2024-01-01" SubmissionCount="1">{significance}</Description>
      </GermlineClassification></RCVClassifications>
    </RCVAccession>"""


VCV = {
    90868: [rcv("RCV9", "Lynch syndrome", "C9", "Pathogenic")],      # de MSH2: no se pide nunca
    101: [rcv("RCV1", "Sickle cell disease", "C1", "Pathogenic"),
          rcv("RCV2", "Malaria, susceptibility to", "C2", "Protective")],
    103: [rcv("RCV3", "beta Thalassemia", "C3", "Pathogenic"),
          rcv("RCV4", "Sickle cell disease", "C1", "Uncertain significance")],
    # 102 y 105: sin RCV en la respuesta (se enlazan sin efecto)
}


def vcv_xml(ids) -> bytes:
    archives = "".join(
        f'<VariationArchive VariationID="{i}" Accession="VCV{i:09d}"><ClassifiedRecord>'
        f'<RCVList>{"".join(VCV[i])}</RCVList></ClassifiedRecord></VariationArchive>'
        for i in ids if i in VCV)
    return f'<?xml version="1.0"?><ClinVarResult-Set>{archives}</ClinVarResult-Set>'.encode()


class FakeNCBI:
    """Responde como el NCBI y apunta cada petición."""

    def __init__(self, docs=None):
        self.docs = CLINVAR_DOCS if docs is None else docs
        self.calls = []

    def __call__(self, url: str) -> bytes:
        parsed = urllib.parse.urlparse(url)
        endpoint = parsed.path.rsplit("/", 1)[-1]
        q = dict(urllib.parse.parse_qsl(parsed.query))
        self.calls.append((endpoint, q))
        assert q["tool"] == "DNAinvent" and q["email"], "faltan tool/email"
        db = q.get("db")
        if db == "gene" and endpoint == "esearch.fcgi":
            return json.dumps({"esearchresult": {"idlist": ["3043"]}}).encode()
        if db == "gene" and endpoint == "esummary.fcgi":
            return json.dumps(GENE_SUMMARY).encode()
        if db == "nuccore" and endpoint == "esearch.fcgi":
            return json.dumps({"esearchresult": {"idlist": ["28302128"]}}).encode()
        if db == "nuccore" and endpoint == "esummary.fcgi":
            return json.dumps({"result": {"uids": ["28302128"],
                                          "28302128": {"accessionversion": "NM_000518.5"}}}).encode()
        if db == "nuccore" and endpoint == "efetch.fcgi":
            assert q["rettype"] == "gbc", "INSDSeq XML se pide con rettype=gbc"
            return insdseq(q["id"])
        if db == "clinvar" and endpoint == "esearch.fcgi":
            assert q["term"] == "HBB[gene]"
            ids = sorted(self.docs)
            start, size = int(q["retstart"]), int(q["retmax"])
            return json.dumps({"esearchresult": {"count": str(len(ids)),
                                                 "idlist": ids[start:start + size]}}).encode()
        if db == "clinvar" and endpoint == "esummary.fcgi":
            uids = q["id"].split(",")
            return json.dumps({"result": {"uids": uids, **{u: self.docs[u] for u in uids}}}).encode()
        if db == "clinvar" and endpoint == "efetch.fcgi":
            assert q["rettype"] == "vcv"
            # Como el ClinVar real: sin is_variationid, los números se leen como otro tipo de ID
            # y llegan OTRAS variantes (aquí, la 90868 de MSH2)
            if q.get("is_variationid") != "true":
                return vcv_xml([90868])
            return vcv_xml(int(i) for i in q["id"].split(","))
        raise AssertionError(f"Petición inesperada: {url}")


@pytest.fixture
def conn():
    connection = open_database(":memory:")
    yield connection
    connection.close()


@pytest.fixture
def fake():
    return FakeNCBI()


@pytest.fixture
def ncbi(fake):
    return NCBIProvider(http=fake, api_key="", sleep=lambda s: None)


@pytest.fixture
def cv(fake):
    return ClinVarProvider(http=fake, api_key="", sleep=lambda s: None)


def count(conn, table):
    return conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]


# ---------------------------------------------------------------- NCBI

def test_parse_location():
    assert parse_location("51..494") == (51, 494)
    assert parse_location("<1..>300") == (1, 300)
    assert parse_location("join(1..10,20..30)") is None


@pytest.mark.parametrize("prefix", ["INSD", "GB"])
def test_parse_both_xml_formats(prefix):
    record = parse_insdseq(insdseq(prefix=prefix), "NM_000518.5")
    assert (record["accession"], record["version"], record["taxon_id"]) == ("NM_000518", 5, 9606)
    assert (record["cds_start"], record["cds_end"]) == (CDS_START, CDS_END)
    assert record["protein"]["residues"] == PROTEIN


def test_unexpected_xml():
    with pytest.raises(ServiceError) as info:
        parse_insdseq(b"<eFetchResult><ERROR>x</ERROR></eFetchResult>", "NM_1.1")
    assert info.value.code == "sequence_not_found"


def test_gene_summary_minus_strand(ncbi):
    summary = ncbi.fetch_gene_summary(3043)
    assert summary["strand"] == -1
    assert (summary["start_pos"], summary["end_pos"]) == (5225464, 5227071)


def test_import_gene(conn, ncbi):
    stages = []
    gene_id = ncbi.import_gene(conn, "HBB", progress=lambda key, data: stages.append(key))
    assert stages == ["database.stage.gene", "database.stage.sequence"]
    gene = queries.find_gene(conn, "hbb")
    assert (gene["id"], gene["ncbi_gene_id"], gene["chromosome"]) == (gene_id, 3043, "11")
    seqs = {s["mol_type"]: s for s in queries.sequences(conn, gene_id)}
    assert (seqs["mrna"]["accession"], seqs["mrna"]["cds_start"]) == ("NM_000518", CDS_START)
    assert seqs["mrna"]["residues"] == MRNA
    assert seqs["protein"]["residues"] == PROTEIN
    ncbi.import_gene(conn, "HBB")                   # dos veces no duplica
    assert (count(conn, "genes"), count(conn, "sequences"), count(conn, "sequence_links")) == (1, 2, 1)


def test_gene_not_found(conn):
    empty = lambda url: json.dumps({"esearchresult": {"idlist": []}}).encode()
    with pytest.raises(ServiceError) as info:
        NCBIProvider(http=empty, api_key="", sleep=lambda s: None).import_gene(conn, "NOEXISTE")
    assert info.value.code == "gene_not_found"


def test_retries_after_429_then_gives_up(fake):
    calls, waits = [], []

    def flaky(url):
        calls.append(url)
        if len(calls) == 1:
            raise urllib.error.HTTPError(url, 429, "Too Many Requests", None, None)
        return fake(url)

    provider = NCBIProvider(http=flaky, api_key="", sleep=waits.append)
    assert provider.search_gene("HBB") == 3043
    assert len(calls) == 2 and 2 in waits

    def down(url):
        raise urllib.error.URLError("sin red")

    with pytest.raises(ServiceError) as info:
        NCBIProvider(http=down, api_key="", sleep=lambda s: None).search_gene("HBB")
    assert info.value.code == "service_unavailable"


def test_api_key_raises_rate(fake):
    assert NCBIProvider(http=fake, api_key="").min_interval > 0.3
    provider = NCBIProvider(http=fake, api_key="ABC", sleep=lambda s: None)
    assert provider.min_interval < 0.2
    provider.search_gene("HBB")
    assert fake.calls[-1][1]["api_key"] == "ABC"


def test_ensure_residues_on_demand(conn, ncbi):
    ncbi.import_gene(conn, "HBB")
    conn.execute("UPDATE sequences SET residues = NULL WHERE mol_type = 'mrna'")
    mrna_id = conn.execute("SELECT id FROM sequences WHERE mol_type = 'mrna'").fetchone()[0]
    assert ncbi.ensure_residues(conn, mrna_id) == MRNA


# ---------------------------------------------------------------- ClinVar: lectura

def test_variation_name():
    v = parse_variation_name("NM_000518.5(HBB):c.20A>T (p.Glu7Val)")
    assert (v["accession"], v["version"], v["hgvs_c"], v["hgvs_p"]) == ("NM_000518", 5, "c.20A>T", "p.Glu7Val")
    assert parse_variation_name("NM_000518.5(HBB):c.93-21G>A")["hgvs_p"] is None
    assert parse_variation_name("GRCh38/hg38 11p15.4(chr11:1-2)x1") is None


def test_real_clinvar_summary():
    lp = parse_summary("5036841", REAL_SAMPLE["5036841"])
    assert (lp["accession"], lp["version"], lp["hgvs_c"], lp["hgvs_p"]) == ("NM_000518", 5, "c.83C>G", "p.Ala28Gly")
    assert lp["clinical_significance"] == "likely pathogenic"
    assert lp["review_status"] == "criteria provided, single submitter"
    assert lp["traits"] == [{"name": "Hb SS disease", "mondo_id": "MONDO:0011382",
                            "medgen_id": "C0002895", "mesh_id": None}]
    lb = parse_summary("5022597", REAL_SAMPLE["5022597"])
    assert (lb["hgvs_p"], lb["traits"]) == ("p.Leu4=", [])     # "Inborn genetic diseases" es genérico
    deletion = parse_summary("5019560", REAL_SAMPLE["5019560"])
    assert (deletion["variant_type"], deletion["ref_allele"]) == ("deletion", None)


def test_somatic_origin():
    v = parse_summary("1", cv_doc("NM_000546.6(TP53):c.524G>A (p.Arg175His)",
                                  "single nucleotide variant", "", somatic="Oncogenic"))
    assert (v["origin"], v["clinical_significance"]) == ("somatic", "oncogenic")


@pytest.mark.parametrize("significance, effect", [
    ("pathogenic/likely pathogenic", "causal"), ("Pathogenic", "causal"),
    ("conflicting classifications of pathogenicity", None), ("likely benign", None),
    ("risk factor", "risk_factor"), ("protective", "protective"), (None, None),
    ("oncogenic", "causal"), ("likely oncogenic", "causal"), ("tier i - strong", None),
    ("established risk allele", "risk_factor"), ("pathogenic, low penetrance", "causal"),
    ("benign/likely benign", None),
])
def test_effects(significance, effect):
    assert classify_effect(significance) == effect


def test_parse_vcv():
    conditions = parse_vcv(vcv_xml([101, 103]))
    assert [(c.rcv, c.name, c.medgen_id, c.significance) for c in conditions[101]] == [
        ("RCV1", "Sickle cell disease", "C1", "pathogenic"),
        ("RCV2", "Malaria, susceptibility to", "C2", "protective")]
    assert conditions[101][0].review_status == "criteria provided, single submitter"
    with pytest.raises(ServiceError):
        parse_vcv(b"esto no es XML")


# ---------------------------------------------------------------- ClinVar: importar

def test_import_variants(conn, cv):
    stages = set()
    stats = cv.import_variants(conn, "HBB", progress=lambda key, data: stages.add(key))
    assert stats.gene_imported
    assert (stats.found, stats.saved, stats.not_on_mrna, stats.other_gene) == (6, 4, 1, 1)
    assert (stats.pathogenic, stats.with_conditions) == (2, 2)
    assert {"database.stage.variants", "database.stage.summaries",
            "database.stage.conditions", "database.stage.saving"} <= stages

    hbs = conn.execute("SELECT * FROM variants WHERE clinvar_id = 101").fetchone()
    assert (hbs["hgvs_p"], hbs["dbsnp_rs"]) == ("p.Glu7Val", "rs334")
    # La variante sobre NM_000518.4 crea solo los metadatos de esa versión
    old = conn.execute("SELECT residues FROM sequences WHERE accession = 'NM_000518' AND version = 4").fetchone()
    assert old["residues"] is None


def test_each_disease_has_its_own_classification(conn, cv):
    """HbS causa la anemia falciforme y protege de la malaria (antes salía como causal)."""
    cv.import_variants(conn, "HBB")
    hbs = conn.execute("SELECT id FROM variants WHERE clinvar_id = 101").fetchone()["id"]
    links = {row["name_en"]: (row["effect"], row["clinical_significance"], row["rcv_accession"])
             for row in queries.variant_diseases(conn, hbs)}
    # El nombre es el del resumen (enlazado por MedGen), así coincide con el resto de variantes
    assert links == {"Hb SS disease": ("causal", "pathogenic", "RCV1"),
                     "Malaria, susceptibility to": ("protective", "protective", "RCV2")}

    gene_id = queries.find_gene(conn, "HBB")["id"]
    table = {row["name_en"]: (row["causal"], row["protective"]) for row in queries.diseases(conn, gene_id)}
    assert table["Hb SS disease"] == (1, 0)          # c.79G>A es incierta para esta enfermedad
    assert table["Malaria, susceptibility to"] == (0, 1)
    causes = {r[0] for r in conn.execute(
        "SELECT d.name_en FROM gene_diseases gd JOIN diseases d ON d.id = gd.disease_id")}
    assert causes == {"Hb SS disease", "beta Thalassemia"}


def test_variant_without_rcv_is_linked_without_effect(conn, cv):
    cv.import_variants(conn, "HBB")
    v = conn.execute("SELECT id FROM variants WHERE clinvar_id = 105").fetchone()["id"]
    assert [(r["name_en"], r["effect"]) for r in queries.variant_diseases(conn, v)] == [("beta Thalassemia", None)]


def test_reimport_fixes_old_links_and_keeps_manual_name(conn, cv, ncbi):
    gene_id = ncbi.import_gene(conn, "HBB")
    mrna_id = queries.reference_mrna(conn, gene_id)["id"]
    conn.execute("INSERT INTO variants (gene_id, sequence_id, hgvs_c, common_name, source_id) "
                 "VALUES (?, ?, 'c.20A>T', 'HbS', 1)", (gene_id, mrna_id))
    # Un enlace antiguo, como los de la primera versión (clasificación global para todo)
    malaria = conn.execute("INSERT INTO diseases (name_en) VALUES ('Malaria, susceptibility to')").lastrowid
    conn.execute("INSERT INTO variant_diseases (variant_id, disease_id, effect, source_id) "
                 "VALUES (1, ?, 'causal', 1)", (malaria,))
    conn.execute("INSERT INTO gene_diseases VALUES (?, ?, 1)", (gene_id, malaria))

    cv.import_variants(conn, "HBB")
    cv.import_variants(conn, "HBB")
    assert count(conn, "variants") == 4
    hbs = conn.execute("SELECT * FROM variants WHERE hgvs_c = 'c.20A>T'").fetchone()
    assert (hbs["clinvar_id"], hbs["common_name"]) == (101, "HbS")
    assert conn.execute("SELECT effect FROM variant_diseases WHERE disease_id = ?", (malaria,)).fetchone()[0] == "protective"
    assert conn.execute("SELECT COUNT(*) FROM gene_diseases WHERE disease_id = ?", (malaria,)).fetchone()[0] == 0


def test_parse_real_vcv():
    """Respuesta real de ClinVar (XML VCV de MSH2 c.2047G>A, noviembre de 2026)."""
    conditions = parse_vcv((DATA / "clinvar_vcv_msh2.xml").read_bytes())
    assert list(conditions) == [90868]
    found = {(c.rcv, c.name, c.medgen_id, c.significance) for c in conditions[90868]}
    assert found == {
        ("RCV000076370", "Lynch syndrome", "C1333990", "pathogenic"),
        ("RCV000132039", "Hereditary cancer-predisposing syndrome", "C0027672", "pathogenic"),
        ("RCV000202225", "not specified", "CN169374", "uncertain significance"),
        ("RCV000524373", "Hereditary nonpolyposis colon cancer", "C0009405", "pathogenic"),
    }
    lynch = next(c for c in conditions[90868] if c.name == "Lynch syndrome")
    assert lynch.review_status == "reviewed by expert panel"


def test_vcv_records_that_were_not_requested_are_ignored(cv):
    stats = clinvar.ImportStats()
    cv.http = lambda url: vcv_xml([101, 90868])        # ClinVar se equivoca y añade otra
    conditions = cv.fetch_conditions([101, 103], stats=stats)
    assert set(conditions) == {101}
    assert (stats.unexpected_records, stats.condition_batches_failed) == (1, 0)


def test_failed_vcv_batch_does_not_stop_import(conn, fake):
    def no_vcv(url):
        if "rettype=vcv" in url:
            raise urllib.error.HTTPError(url, 400, "Bad Request", None, None)
        return fake(url)

    stats = ClinVarProvider(http=no_vcv, api_key="", sleep=lambda s: None).import_variants(conn, "HBB")
    assert (stats.saved, stats.with_conditions, stats.condition_batches_failed) == (4, 0, 1)
    assert count(conn, "gene_diseases") == 0          # sin RCV no se afirma que cause nada
    assert conn.execute("SELECT COUNT(*) FROM variant_diseases WHERE effect IS NOT NULL").fetchone()[0] == 0


def test_only_pathogenic(conn, cv):
    stats = cv.import_variants(conn, "HBB", only_pathogenic=True)
    assert (stats.saved, stats.filtered) == (2, 2)


def test_paging_and_max(cv, monkeypatch):
    monkeypatch.setattr(clinvar, "SEARCH_PAGE", 2)
    monkeypatch.setattr(clinvar, "SUMMARY_BATCH", 2)
    assert len(cv.search_variant_ids("HBB")) == 6
    assert len(cv.search_variant_ids("HBB", max_variants=3)) == 3
    assert len(cv.fetch_summaries(["101", "102", "103"])) == 3


# ---------------------------------------------------------------- variantes y consultas

def imported(conn, cv):
    cv.import_variants(conn, "HBB")
    return queries.find_gene(conn, "HBB")["id"]


def test_sickle_cell_effect_and_side_by_side(conn, cv):
    gene_id = imported(conn, cv)
    [hbs] = queries.find_variants(conn, gene_id, "rs334")
    effect, residues, cds_start = locate(conn, hbs["id"])
    assert (effect.position, effect.codon_number) == (CDS_START + 19, 7)
    assert (effect.ref_codon, effect.alt_codon, effect.kind, effect.hgvs_p) == ("GAG", "GTG", "missense", "p.Glu7Val")
    view = side_by_side(effect, residues, cds_start, context=2)
    assert view.first_codon == 5
    assert view.healthy == ["ACT", "CCT", "GAG", "GAG", "AAG"]
    assert view.variant == ["ACT", "CCT", "GTG", "GAG", "AAG"]
    assert (view.healthy_protein, view.variant_protein) == ("TPEEK", "TPVEK")
    assert (view.changed_codon, view.changed_base) == (2, 1)


def test_variant_checks(conn, cv):
    gene_id = imported(conn, cv)
    source = conn.execute("SELECT id FROM sources LIMIT 1").fetchone()[0]
    mrna = queries.reference_mrna(conn, gene_id)["id"]

    def add(hgvs):
        return conn.execute("INSERT INTO variants (gene_id, sequence_id, hgvs_c, source_id) "
                            "VALUES (?, ?, ?, ?)", (gene_id, mrna, hgvs, source)).lastrowid

    for hgvs, code in [("c.21A>T", "variant_mismatch"), ("c.999A>T", "variant_outside_cds"),
                       ("c.93-21G>A", "variant_not_supported")]:
        with pytest.raises(InvalidParameterError) as info:
            locate(conn, add(hgvs))
        assert info.value.code == code
    effect, _, _ = locate(conn, add("c.21G>A"))           # GAG -> GAA: también ácido glutámico
    assert (effect.kind, effect.hgvs_p) == ("synonymous", "p.Glu7=")
    effect, _, _ = locate(conn, add("c.19G>T"))           # GAG -> TAG: parada
    assert (effect.kind, effect.hgvs_p) == ("nonsense", "p.Glu7Ter")


def test_find_variants(conn, cv):
    gene_id = imported(conn, cv)
    for text in ["c.20A>T", "p.Glu7Val", "Glu7Val", "rs334", "RS334"]:
        assert [v["clinvar_id"] for v in queries.find_variants(conn, gene_id, text)] == [101], text
    assert len(queries.find_variants(conn, gene_id, "c.20")) == 2         # c.20A>T y c.20A>G
    assert [v["clinvar_id"] for v in queries.find_variants(conn, gene_id, causal_only=True)] == [101, 103]


def test_significance_counts_order(conn, cv):
    gene_id = imported(conn, cv)
    assert [s for s, _ in queries.significance_counts(conn, gene_id)] == [
        "pathogenic/likely pathogenic", "conflicting classifications of pathogenicity",
        "likely benign", "pathogenic; other"]


def test_search(conn, cv):
    imported(conn, cv)
    assert [r["kind"] for r in queries.search(conn, "hemoglobin")] == ["gene"]
    labels = [r["label"] for r in queries.search(conn, "thalass")]
    assert labels == ["beta Thalassemia"]
    assert queries.search(conn, "  ") == []
    assert queries.search(conn, '"malaria*') != []            # comillas y asteriscos no rompen FTS


def test_migration_from_version_1(tmp_path):
    """Una base de datos creada con la primera versión (sin columnas por enfermedad)."""
    path = tmp_path / "old.db"
    old_schema = SCHEMA_PATH.read_text(encoding="utf-8")
    old_schema = (old_schema.replace("    review_status   TEXT,                      -- revisión de la clasificación global (ClinVar)\n", "")
                  .replace("    clinical_significance TEXT,                -- 'pathogenic', 'protective'... para esta enfermedad\n", "")
                  .replace("    review_status   TEXT,                      -- 'criteria provided, single submitter'...\n", "")
                  .replace("    rcv_accession   TEXT,                      -- 'RCV000016574'\n", "")
                  .replace("PRAGMA user_version = 4;", ""))
    raw = sqlite3.connect(path)
    raw.executescript(old_schema)
    raw.close()
    connection = open_database(path)
    columns = {r[1] for r in connection.execute("PRAGMA table_info(variant_diseases)")}
    assert {"clinical_significance", "review_status", "rcv_accession"} <= columns
    assert connection.execute("PRAGMA user_version").fetchone()[0] == 4
    connection.close()
    open_database(path).close()                       # abrirla otra vez no falla


def test_console_dynamic_keys_exist_in_every_language():
    """Las claves que la consola construye con f-strings (el test de i18n no las ve)."""
    from i18n import Translator, available_languages
    from ui import console

    keys = [f"database.significance.{k}" for k in set(console.SIGNIFICANCE_KEYS.values()) | {"none"}]
    keys += [f"database.effects.{e}" for e in ("causal", "protective", "risk_factor", "none")]
    keys += [f"database.variant_kinds.{k}" for k in ("missense", "synonymous", "nonsense", "stop_lost")]
    keys += [f"database.mol_types.{m}" for m in ("mrna", "protein", "genomic", "cds", "ncrna")]
    keys += [f"units.{u}.{n}" for u in ("variant", "variant_causal", "variant_protective",
                                        "variant_risk_factor", "atom") for n in ("one", "other")]
    keys += [f"database.phase.{k}" for k in ("approved", "phase3", "phase2", "phase1", "early", "unknown")]
    keys += [f"database.{group}.{s}" for group in ("status", "status_short")
             for s in ("approved", "withdrawn", "trials", "none")]
    keys += [f"database.ligand_kinds.{k}" for k in ("drug", "cofactor", "additive")]
    keys += [f"database.contact_kinds.{k}" for k in ("covalent", "hbond", "hydrophobic", "contact")]
    keys += [f"database.actions.{k}" for k in set(console.ACTION_KEYS.values())]
    for language in available_languages():
        tr = Translator(language)
        tr._fallback = {}                     # sin el español de reserva: que falte se note
        missing = [key for key in keys if tr(key) == key]
        assert not missing, (language, missing)
