"""NCBIProvider.gene_overview: la ficha básica de un gen que no es de serie, sin red."""

import json
import urllib.parse

from dnainvent.database.ncbi import NCBIProvider


def fake_ncbi(asked):
    def http(url):
        params = dict(urllib.parse.parse_qsl(urllib.parse.urlparse(url).query))
        asked.append((url.split("/")[-1].split("?")[0], params.get("db"), params.get("term") or params.get("id")))
        if params["db"] == "gene" and "esearch" in url:
            ids = ["4625"] if params["term"].startswith(("MYH7[", "CMH1[")) and "[sym]" in params["term"] or params["term"].startswith("MYH7[pref]") else []
            return json.dumps({"esearchresult": {"idlist": ids}}).encode()
        if params["db"] == "gene" and "esummary" in url:
            if "," in params["id"]:
                return json.dumps({"result": {"uids": ["4625"], "4625": {"name": "MYH7"}}}).encode()
            return json.dumps({"result": {"4625": {
                "name": "MYH7", "description": "myosin heavy chain 7", "chromosome": "14",
                "maplocation": "14q11.2", "otheraliases": "CMD1S, CMH1, MPD1", "mim": ["160760"],
                "summary": "Muscle myosin.", "organism": {"taxid": 9606}}}}).encode()
        if params["db"] == "clinvar":
            return json.dumps({"esearchresult": {"count": "4321", "idlist": []}}).encode()
        raise AssertionError(url)
    return http


def provider(asked):
    return NCBIProvider(http=fake_ncbi(asked), api_key="", sleep=lambda s: None)


def test_gene_overview():
    asked = []
    gene = provider(asked).gene_overview("MYH7")
    assert gene == {"ncbi_gene_id": 4625, "symbol": "MYH7", "name": "myosin heavy chain 7",
                    "chromosome": "14", "map_location": "14q11.2", "aliases": ["CMD1S", "CMH1", "MPD1"],
                    "summary": "Muscle myosin.", "mim": ["160760"], "taxon_id": 9606,
                    "clinvar_records": 4321}
    assert ("esearch.fcgi", "clinvar", "MYH7[gene]") in asked


def test_gene_overview_by_alias_gives_the_official_symbol():
    gene = provider([]).gene_overview("CMH1")
    assert gene["symbol"] == "MYH7"


def test_gene_overview_of_an_unknown_gene():
    asked = []
    assert provider(asked).gene_overview("NOEXISTE") is None
    assert all(db != "clinvar" for _, db, _ in asked)
