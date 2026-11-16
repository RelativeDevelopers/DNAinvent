"""Tests del catálogo de genes de serie y de Orphanet (sin red).

Los XML de Orphanet son de prueba, con la estructura que describe «Orphadata Science Data
Description» (los reales ocupan varios MB): el código 232 es la anemia falciforme y el 848 la
beta talasemia, como en Orphanet, pero los demás datos son inventados para el test.
"""

import pytest

from dnainvent import ServiceError
from dnainvent.database import open_database
from dnainvent.database.catalog import COMMON_NAMES, CatalogBuilder, catalog_genes
from dnainvent.database.chembl import DrugStats
from dnainvent.database.clinvar import ClinVarProvider
from dnainvent.database.ncbi import NCBIProvider
from dnainvent.database.orphanet import (OrphanetProvider, Prevalence, best_prevalence,
                                         parse_natural_history, parse_prevalence)
from dnainvent.protein.cache import Cache
from test_database import CLINVAR_DOCS, FakeNCBI, cv_doc

AGES = """<?xml version="1.0" encoding="UTF-8"?>
<JDBOR date="2026-06-23 07:57:58" version="1.3.42" copyright="Orphanet (c) 2026">
  <Availability><Licence><FullName lang="en">Creative Commons Attribution 4.0 International</FullName>
    <ShortIdentifier>CC-BY-4.0</ShortIdentifier></Licence></Availability>
  <DisorderList count="3">
    <Disorder id="1"><OrphaCode>232</OrphaCode><Name lang="en">Sickle cell anemia</Name>
      <AverageAgeOfOnsetList count="2">
        <AverageAgeOfOnset id="1"><Name lang="en">Infancy</Name></AverageAgeOfOnset>
        <AverageAgeOfOnset id="2"><Name lang="en">Childhood</Name></AverageAgeOfOnset>
      </AverageAgeOfOnsetList>
      <TypeOfInheritanceList count="2">
        <TypeOfInheritance id="1"><Name lang="en">Autosomal recessive</Name></TypeOfInheritance>
        <TypeOfInheritance id="2"><Name lang="en">Autosomal recessive</Name></TypeOfInheritance>
      </TypeOfInheritanceList>
    </Disorder>
    <Disorder id="2"><OrphaCode>848</OrphaCode><Name lang="en">beta Thalassemia</Name>
      <TypeOfInheritanceList count="1">
        <TypeOfInheritance id="1"><Name lang="en">Autosomal recessive</Name></TypeOfInheritance>
      </TypeOfInheritanceList>
    </Disorder>
    <Disorder id="3"><OrphaCode>999</OrphaCode><Name lang="en">Unknown thing</Name>
      <AverageAgeOfOnsetList count="1">
        <AverageAgeOfOnset id="1"><Name lang="en">No data available</Name></AverageAgeOfOnset>
      </AverageAgeOfOnsetList>
    </Disorder>
  </DisorderList>
</JDBOR>""".encode()


def prevalence(kind, qualification, klass, geo, status, value="0.0"):
    return f"""<Prevalence id="1"><Source>Test</Source>
      <PrevalenceType id="1"><Name lang="en">{kind}</Name></PrevalenceType>
      <PrevalenceQualification id="1"><Name lang="en">{qualification}</Name></PrevalenceQualification>
      <PrevalenceClass id="1"><Name lang="en">{klass}</Name></PrevalenceClass>
      <ValMoy>{value}</ValMoy>
      <PrevalenceGeographic id="1"><Name lang="en">{geo}</Name></PrevalenceGeographic>
      <PrevalenceValidationStatus id="1"><Name lang="en">{status}</Name></PrevalenceValidationStatus>
    </Prevalence>"""


PREV = f"""<?xml version="1.0" encoding="UTF-8"?>
<JDBOR date="2026-06-23 08:00:00" version="1.3.42">
  <DisorderList count="2">
    <Disorder id="1"><OrphaCode>232</OrphaCode><Name lang="en">Sickle cell anemia</Name>
      <PrevalenceList count="4">
        {prevalence("Point prevalence", "Value and class", "1-5 / 10 000", "Spain", "Validated")}
        {prevalence("Point prevalence", "Value and class", "6-9 / 10 000", "Worldwide", "Not yet validated")}
        {prevalence("Birth prevalence", "Value and class", "1-5 / 10 000", "Europe", "Validated")}
        {prevalence("Cases/families", "Case", "", "Worldwide", "Validated", value="12")}
      </PrevalenceList>
    </Disorder>
    <Disorder id="2"><OrphaCode>848</OrphaCode><Name lang="en">beta Thalassemia</Name>
      <PrevalenceList count="1">
        {prevalence("Point prevalence", "Class only", "Unknown", "Worldwide", "Validated")}
      </PrevalenceList>
    </Disorder>
  </DisorderList>
</JDBOR>""".encode()


def orphanet_http(calls=None, broken_first=False):
    def http(url):
        if calls is not None:
            calls.append(url)
        if broken_first and "orphadata.com" in url:
            return b"version https://git-lfs.github.com/spec/v1\noid sha256:abc\n"   # no es el XML
        return AGES if "ages" in url else PREV
    return http


# ---------------------------------------------------------------- Orphanet

def test_parse_natural_history():
    history, date = parse_natural_history(AGES)
    assert date == "2026-06-23 07:57:58"
    assert history["232"].inheritance == ["Autosomal recessive"]          # sin repetir
    assert history["232"].onset == ["Infancy", "Childhood"]
    assert history["999"].onset == []                                     # «No data available»


def test_best_prevalence_prefers_validated_and_wide_areas():
    found, _ = parse_prevalence(PREV)
    name, items = found["232"]
    best = best_prevalence(items)
    # Validada > sin validar; Europa > un país; puntual > de nacimiento
    assert (best.type, best.prevalence_class, best.geographic) == ("Birth prevalence", "1-5 / 10 000", "Europe")
    assert best_prevalence(found["848"][1]) is None                       # clase «Unknown»
    assert best_prevalence([]) is None


def test_point_prevalence_beats_birth_prevalence_in_the_same_area():
    items = [Prevalence("Birth prevalence", "Class only", "1-9 / 100 000", "Worldwide", "Validated", "", ""),
             Prevalence("Point prevalence", "Class only", "1-5 / 10 000", "Worldwide", "Validated", "", "")]
    assert best_prevalence(items).type == "Point prevalence"


def test_fetch_falls_back_to_github_and_caches(tmp_path):
    calls = []
    provider = OrphanetProvider(http=orphanet_http(calls, broken_first=True), cache=Cache(tmp_path))
    assert provider.fetch("ages") == AGES
    assert ["orphadata.com" in calls[0], "media.githubusercontent.com" in calls[1]] == [True, True]
    provider.fetch("ages")
    assert len(calls) == 2                                                # segunda vez, de la caché


def test_fetch_from_local_folder(tmp_path):
    (tmp_path / "en_product9_ages.xml").write_bytes(AGES)
    provider = OrphanetProvider(http=lambda url: pytest.fail("no debería descargar"),
                                cache=Cache(tmp_path / "cache"), local_dir=tmp_path)
    assert provider.fetch("ages") == AGES


def test_fetch_fails_cleanly(tmp_path):
    def down(url):
        raise OSError("sin red")
    with pytest.raises(ServiceError) as info:
        OrphanetProvider(http=down, cache=Cache(tmp_path)).fetch("prev")
    assert info.value.code == "service_unavailable"


# ---------------------------------------------------------------- catálogo

class FakeChEMBL:
    def __init__(self, fail_diseases=False):
        self.fail_diseases = fail_diseases
        self.gene_calls, self.disease_calls = [], []

    def import_gene_drugs(self, conn, gene_id, progress=None):
        self.gene_calls.append(gene_id)
        return DrugStats(accession="P68871", drugs=2)

    def import_disease_drugs(self, conn, disease_id, progress=None):
        self.disease_calls.append(disease_id)
        if self.fail_diseases:
            raise ServiceError("service_unavailable", service="ChEMBL")
        return DrugStats()


@pytest.fixture
def fake():
    docs = dict(CLINVAR_DOCS)
    # HbS, con los códigos de Orphanet que da ClinVar para la anemia falciforme
    docs["101"] = cv_doc("NM_000518.5(HBB):c.20A>T (p.Glu7Val)", "single nucleotide variant",
                         "Pathogenic", [("Hb SS disease", [{"db_source": "MedGen", "db_id": "C1"},
                                                           {"db_source": "Orphanet", "db_id": "232"}])],
                         rs="334")
    return FakeNCBI(docs)


def make_builder(conn, fake, tmp_path, chembl=None):
    ncbi = NCBIProvider(http=fake, api_key="", sleep=lambda s: None)
    return CatalogBuilder(conn, ncbi=ncbi,
                          clinvar=ClinVarProvider(http=fake, api_key="", sleep=lambda s: None),
                          chembl=chembl or FakeChEMBL(),
                          orphanet=OrphanetProvider(http=orphanet_http(), cache=Cache(tmp_path)))


def test_build_catalog(tmp_path, fake):
    conn = open_database(tmp_path / "catalog.db")
    chembl = FakeChEMBL()
    report = make_builder(conn, fake, tmp_path, chembl).run(["HBB"])
    assert not report.failed, report.failed
    assert [(s.step, s.status) for s in report.steps] == [
        ("gene", "done"), ("variants", "done"), ("sequences", "done"), ("gene_drugs", "done"),
        ("disease_drugs", "done"), ("orphanet", "done"), ("common_names", "done")]
    assert catalog_genes(conn) == ["HBB"]

    # Nombre conocido; los de otros genes no están en este catálogo y no se cuentan
    assert conn.execute("SELECT common_name FROM variants WHERE hgvs_c = 'c.20A>T'").fetchone()[0] == "HbS"
    assert ("HBB", "HbC") in report.missing_names
    assert all(gene == "HBB" for gene, _ in report.missing_names)

    # El ARNm antiguo (NM_000518.4, de la variante 103) se ha descargado para comparar sin conexión
    pending = conn.execute("SELECT COUNT(*) FROM sequences s JOIN variants v ON v.sequence_id = s.id "
                           "WHERE s.residues IS NULL").fetchone()[0]
    assert pending == 0

    # Orphanet: por código (anemia falciforme) y por nombre (beta talasemia, sin código en ClinVar)
    sickle = conn.execute("SELECT * FROM diseases WHERE name_en = 'Hb SS disease'").fetchone()
    assert (sickle["inheritance"], sickle["onset"]) == ("Autosomal recessive", "Infancy|Childhood")
    assert (sickle["prevalence_class"], sickle["prevalence_geo"]) == ("1-5 / 10 000", "Europe")
    thal = conn.execute("SELECT * FROM diseases WHERE name_en = 'beta Thalassemia'").fetchone()
    assert (thal["orpha_codes"], thal["inheritance"], thal["prevalence_class"]) == ("848", "Autosomal recessive", None)

    info = dict(conn.execute("SELECT key, value FROM catalog_info").fetchall())
    assert info["genes"] == "HBB" and info["orphadata_ages_date"] == "2026-06-23 07:57:58"
    assert chembl.disease_calls                         # ha buscado fármacos de sus enfermedades
    conn.close()


def test_build_resumes_where_it_stopped(tmp_path, fake):
    path = tmp_path / "catalog.db"
    conn = open_database(path)
    report = make_builder(conn, fake, tmp_path, FakeChEMBL(fail_diseases=True)).run(["HBB"])
    assert [(s.step, s.detail) for s in report.failed] == [
        ("disease_drugs", "2 de 2 enfermedades sin respuesta de ChEMBL")]
    assert catalog_genes(conn) == []                    # le falta un paso
    conn.close()

    # Segunda vez: lo hecho no se repite (ninguna consulta al NCBI) y lo fallido sí
    fake.calls.clear()
    conn = open_database(path)
    chembl = FakeChEMBL()
    report = make_builder(conn, fake, tmp_path, chembl).run(["HBB"])
    statuses = {s.step: s.status for s in report.steps}
    assert statuses == {"gene": "skipped", "variants": "skipped", "sequences": "skipped",
                        "gene_drugs": "skipped", "disease_drugs": "done",
                        "orphanet": "done", "common_names": "done"}
    assert fake.calls == [] and chembl.gene_calls == []
    assert catalog_genes(conn) == ["HBB"]
    conn.close()


def test_unknown_gene_stops_only_that_gene(tmp_path, fake):
    conn = open_database(tmp_path / "catalog.db")

    def no_gene(url):
        if "db=gene" in url and "NOEXISTE" in url:
            return b'{"esearchresult": {"idlist": []}}'
        return fake(url)

    ncbi = NCBIProvider(http=no_gene, api_key="", sleep=lambda s: None)
    builder = CatalogBuilder(conn, ncbi=ncbi,
                             clinvar=ClinVarProvider(http=no_gene, api_key="", sleep=lambda s: None),
                             chembl=FakeChEMBL(),
                             orphanet=OrphanetProvider(http=orphanet_http(), cache=Cache(tmp_path)))
    report = builder.run(["NOEXISTE", "HBB"])
    noexiste = [s for s in report.steps if s.symbol == "NOEXISTE"]
    assert [(s.step, s.status) for s in noexiste] == [("gene", "failed")]
    assert "gene_not_found" in noexiste[0].detail
    assert catalog_genes(conn) == ["HBB"]
    conn.close()


def test_common_names_are_unique_and_well_formed():
    keys = [(gene, c) for gene, c, _, _ in COMMON_NAMES]
    assert len(keys) == len(set(keys))
    assert all(c.startswith("c.") and p.startswith("p.") for _, c, p, _ in COMMON_NAMES)


# ---------------------------------------------------------------- fallos de la primera ejecución real

def gene_search_http(pref_ids, sym_ids, names):
    """NCBI Gene de prueba: [pref] y [sym] devuelven listas distintas; esummary, los símbolos."""
    import json
    from urllib.parse import parse_qs, urlparse

    def http(url):
        q = {k: v[0] for k, v in parse_qs(urlparse(url).query).items()}
        if "esearch" in url:
            ids = pref_ids if "[pref]" in q["term"] else sym_ids
            return json.dumps({"esearchresult": {"idlist": ids}}).encode()
        uids = q["id"].split(",")
        return json.dumps({"result": {"uids": uids, **{u: {"name": names[u]} for u in uids}}}).encode()
    return http


def test_search_gene_uses_the_official_symbol():
    http = gene_search_http(pref_ids=["3064"], sym_ids=["6532", "3064"],
                            names={"6532": "SLC6A4", "3064": "HTT"})
    assert NCBIProvider(http=http, api_key="", sleep=lambda s: None).search_gene("HTT") == 3064


def test_search_gene_with_aliases_picks_the_exact_symbol():
    """Si [pref] no da nada, de los resultados de [sym] se coge el del símbolo exacto:
    «HTT» también es un alias de SLC6A4 (5-HTT), que llegaba el primero."""
    http = gene_search_http(pref_ids=[], sym_ids=["6532", "3064"],
                            names={"6532": "SLC6A4", "3064": "HTT"})
    assert NCBIProvider(http=http, api_key="", sleep=lambda s: None).search_gene("HTT") == 3064


def test_failed_vcv_batch_is_retried_in_small_pieces(fake):
    """BRCA1: 84 lotes de 25 fallaron. Un lote que falla se repite en trozos de 5."""
    from test_database import vcv_xml
    import urllib.error
    from urllib.parse import parse_qs, urlparse

    calls = []

    def http(url):
        q = {k: v[0] for k, v in parse_qs(urlparse(url).query).items()}
        if q.get("rettype") == "vcv":
            ids = q["id"].split(",")
            calls.append(len(ids))
            if len(ids) > 5:
                raise urllib.error.URLError("respuesta demasiado grande")
            return vcv_xml(ids)
        return fake(url)

    provider = ClinVarProvider(http=http, api_key="", sleep=lambda s: None)
    from dnainvent.database.clinvar import ImportStats
    stats = ImportStats()
    ids = list(range(1, 13))
    provider.fetch_conditions(ids, stats=stats)
    assert stats.condition_batches_failed == 0
    assert [n for n in calls if n <= 5] == [5, 5, 2]          # los tres trozos del lote de 12


def test_wrong_gene_from_ncbi_is_not_kept(tmp_path, fake):
    """Si el NCBI devuelve otro gen (antes: SLC6A4 al pedir HTT), el paso falla y no se guarda."""
    conn = open_database(tmp_path / "catalog.db")
    report = make_builder(conn, fake, tmp_path).run(["HTT"])    # FakeNCBI siempre devuelve HBB
    htt = [s for s in report.steps if s.symbol == "HTT"]
    assert [(s.step, s.status) for s in htt] == [("gene", "failed")]
    assert "HBB en lugar de HTT" in htt[0].detail
    assert conn.execute("SELECT COUNT(*) FROM genes").fetchone()[0] == 0
    conn.close()


def test_clean_removes_foreign_genes_and_forgets_their_steps(tmp_path, fake):
    """Repara un catálogo como el de la primera ejecución: SLC6A4 guardado y los pasos de HTT
    apuntados como hechos."""
    conn = open_database(tmp_path / "catalog.db")
    builder = make_builder(conn, fake, tmp_path)
    builder.run(["HBB"])
    conn.execute("UPDATE genes SET symbol = 'SLC6A4' WHERE symbol = 'HBB'")
    conn.execute("UPDATE catalog_steps SET symbol = 'HTT' WHERE symbol = 'HBB'")
    conn.commit()

    assert builder.clean() == ["SLC6A4"]
    for table in ("genes", "variants", "variant_diseases", "sequences", "diseases", "catalog_steps"):
        expected = 2 if table == "catalog_steps" else 0          # quedan los dos pasos de '*'
        assert conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0] == expected, table
    assert conn.execute("SELECT COUNT(*) FROM search_index WHERE kind IN ('gene', 'disease')").fetchone()[0] == 0
    conn.close()


def test_unexpected_error_fails_only_its_step(tmp_path, fake):
    class BrokenChEMBL(FakeChEMBL):
        def import_gene_drugs(self, conn, gene_id, progress=None):
            raise TypeError("'NoneType' object is not subscriptable")

    conn = open_database(tmp_path / "catalog.db")
    report = make_builder(conn, fake, tmp_path, BrokenChEMBL()).run(["HBB"])
    statuses = {s.step: s.status for s in report.steps}
    assert statuses["gene_drugs"] == "failed" and statuses["disease_drugs"] == "done"
    assert "error interno: TypeError" in report.failed[0].detail
    conn.close()
