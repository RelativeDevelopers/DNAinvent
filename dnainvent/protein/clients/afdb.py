"""AlphaFold DB (EMBL-EBI y Google DeepMind): estructuras predichas por AlphaFold.

    db = AlphaFoldDB()
    entry, structure = db.structure("P68871")       # hemoglobina beta humana
    entry.model_id    -> "AF-P68871-F1"
    entry.version     -> 6

1. Ficha de la entrada:  https://alphafold.ebi.ac.uk/api/prediction/P68871
   Devuelve una lista JSON (puede incluir complejos de varias cadenas). Se elige el
   modelo de la cadena sola, "AF-<accesión>-F1".
2. Estructura en PDB: la URL "pdbUrl" de la ficha.

En la v6 (2025) la API renombró campos (entryId -> modelEntityId, uniprotSequence ->
sequence...) y los antiguos dejaron de existir el 25-06-2026. Se aceptan los dos nombres.

Los datos de AlphaFold DB tienen licencia CC-BY 4.0: hay que citarlos al mostrarlos
(Jumper et al. 2021, Nature; Varadi et al. 2024, Nucleic Acids Research).
"""

import json
import re
import urllib.error
from dataclasses import dataclass
from typing import Callable

from ...errors import InvalidParameterError, ServiceError
from ..cache import Cache
from ..pdb import Structure, parse_pdb
from . import http_get

SERVICE = "AlphaFold DB"
API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
FILES_URL = "https://alphafold.ebi.ac.uk/files/{model_id}-model_v{version}.pdb"
ENTRY_MAX_AGE = 30 * 24 * 3600          # la ficha se renueva cada 30 días (versiones nuevas)

# Formato oficial de las accesiones de UniProt (P68871, Q9Y6K9, A0A023GPI8...)
ACCESSION = re.compile(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}")


@dataclass(frozen=True)
class AlphaFoldEntry:
    accession: str
    model_id: str                 # "AF-P68871-F1"
    version: int | None           # versión del modelo (6 en 2025-2026)
    sequence: str                 # secuencia de UniProt que se plegó
    pdb_url: str
    mean_plddt: float | None      # confianza media (0-100)

    @classmethod
    def from_json(cls, accession: str, item: dict) -> "AlphaFoldEntry":
        model_id = item.get("modelEntityId") or item.get("entryId") or f"AF-{accession}-F1"
        version = item.get("latestVersion")
        pdb_url = item.get("pdbUrl") or FILES_URL.format(model_id=model_id, version=version)
        return cls(
            accession=accession, model_id=model_id,
            version=int(version) if version is not None else None,
            sequence=item.get("sequence") or item.get("uniprotSequence") or "",
            pdb_url=pdb_url, mean_plddt=item.get("globalMetricValue"),
        )


def validate_accession(accession: str) -> str:
    accession = accession.strip().upper()
    if not ACCESSION.fullmatch(accession):
        raise InvalidParameterError("invalid_accession", accession=accession)
    return accession


class AlphaFoldDB:
    def __init__(self, cache: Cache | None = None, http: Callable[[str], bytes] = http_get):
        self.cache = cache if cache is not None else Cache()
        self.http = http

    def entry(self, accession: str) -> AlphaFoldEntry:
        accession = validate_accession(accession)
        raw = self.cache.get("afdb", f"{accession}.json", max_age=ENTRY_MAX_AGE)
        if raw is None:
            raw = self._download(API_URL.format(accession=accession), accession)
        try:
            items = json.loads(raw)
        except ValueError:
            raise ServiceError("service_unavailable", service=SERVICE) from None
        item = _choose(items, accession)
        if item is None:
            raise ServiceError("protein_not_found", accession=accession)
        self.cache.put("afdb", f"{accession}.json", raw)
        return AlphaFoldEntry.from_json(accession, item)

    def structure_text(self, accession: str) -> tuple[AlphaFoldEntry, str]:
        entry = self.entry(accession)
        key = f"{entry.model_id}-v{entry.version}.pdb"
        raw = self.cache.get("afdb", key)
        if raw is None:
            raw = self._download(entry.pdb_url, accession)
            self.cache.put("afdb", key, raw)
        return entry, raw.decode("utf-8", errors="replace")

    def structure(self, accession: str) -> tuple[AlphaFoldEntry, Structure]:
        entry, text = self.structure_text(accession)
        return entry, parse_pdb(text)

    def _download(self, url: str, accession: str) -> bytes:
        try:
            return self.http(url)
        except urllib.error.HTTPError as error:
            if error.code == 404:
                raise ServiceError("protein_not_found", accession=accession) from None
            raise ServiceError("service_unavailable", service=SERVICE) from None
        except (urllib.error.URLError, OSError):       # sin red, DNS, tiempo agotado...
            raise ServiceError("service_unavailable", service=SERVICE) from None


def _choose(items, accession: str) -> dict | None:
    """De la lista de la API, el modelo de la cadena sola (F1) de esta accesión."""
    if not isinstance(items, list):
        return None
    own = [item for item in items if isinstance(item, dict)
           and str(item.get("modelEntityId") or item.get("entryId") or "")
           .startswith(f"AF-{accession}-F")]
    preferred = [item for item in own
                 if (item.get("modelEntityId") or item.get("entryId")) == f"AF-{accession}-F1"]
    return (preferred or own or [None])[0]
