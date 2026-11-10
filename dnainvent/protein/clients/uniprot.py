"""UniProt: de un gen a su proteína (la accesión que usan ChEMBL, AlphaFold DB y el PDB).

    UniProt().accession("HBB")           -> "P68871"
    UniProt().accession("Hbb-b1", 10090) -> la del ratón

Busca entradas revisadas (Swiss-Prot) con ese nombre de gen exacto en ese organismo:
https://rest.uniprot.org/uniprotkb/search?query=gene_exact:HBB AND organism_id:9606 AND
reviewed:true&fields=accession,gene_primary&format=json

Datos de UniProt: CC-BY 4.0.
"""

import json
import urllib.error
import urllib.parse
from typing import Callable

from ...errors import ServiceError
from ..cache import Cache
from . import http_get

SERVICE = "UniProt"
SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search?"
MAX_AGE = 30 * 24 * 3600


class UniProt:
    def __init__(self, cache: Cache | None = None, http: Callable[[str], bytes] = http_get):
        self.cache = cache if cache is not None else Cache()
        self.http = http

    def accession(self, symbol: str, taxon_id: int = 9606) -> str | None:
        """Accesión de la proteína revisada de ese gen, o None si no hay ninguna."""
        key = f"{symbol.upper()}_{taxon_id}.json"
        raw = self.cache.get("uniprot", key, max_age=MAX_AGE)
        if raw is None:
            query = f"gene_exact:{symbol} AND organism_id:{taxon_id} AND reviewed:true"
            url = SEARCH_URL + urllib.parse.urlencode(
                {"query": query, "fields": "accession,gene_primary", "format": "json", "size": 5})
            try:
                raw = self.http(url)
            except (urllib.error.URLError, OSError):
                raise ServiceError("service_unavailable", service=SERVICE) from None
            self.cache.put("uniprot", key, raw)
        try:
            results = json.loads(raw).get("results", [])
        except ValueError:
            raise ServiceError("service_unavailable", service=SERVICE) from None
        # Primero la que tiene ese gen como nombre principal (no como sinónimo)
        for item in results:
            genes = item.get("genes") or [{}]
            primary = (genes[0].get("geneName") or {}).get("value", "")
            if primary.upper() == symbol.upper():
                return item.get("primaryAccession")
        return results[0].get("primaryAccession") if results else None
