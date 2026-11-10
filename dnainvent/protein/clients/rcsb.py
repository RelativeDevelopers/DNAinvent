"""RCSB PDB: estructuras medidas en el laboratorio (cristalografía, crio-EM, RMN).

    pdb = RCSB()
    structure = pdb.structure("5E83")    # hemoglobina S (anemia falciforme) con voxelotor

A diferencia de AlphaFold DB, aquí hay estructuras con el fármaco DENTRO de la proteína, tal
como se midió. Se descargan en formato PDB de https://files.rcsb.org/download/<ID>.pdb y se
guardan en la caché sin caducidad (una entrada publicada no cambia; si se corrige, cambia de
versión pero mantiene el ID y basta con borrar la caché).

Las estructuras muy grandes (ribosomas, virus) solo existen en mmCIF: entonces se avisa con
"structure_not_found".

Datos del PDB: dominio público (CC0). Al mostrarlos se cita el ID de la entrada.
"""

import re
import urllib.error
from typing import Callable

from ...errors import InvalidParameterError, ServiceError
from ..cache import Cache
from ..pdb import Structure, parse_pdb
from . import http_get

SERVICE = "RCSB PDB"
FILE_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
PDB_ID = re.compile(r"[1-9][A-Z0-9]{3}")


def validate_pdb_id(pdb_id: str) -> str:
    pdb_id = pdb_id.strip().upper()
    if not PDB_ID.fullmatch(pdb_id):
        raise InvalidParameterError("invalid_pdb_id", pdb_id=pdb_id)
    return pdb_id


class RCSB:
    def __init__(self, cache: Cache | None = None, http: Callable[[str], bytes] = http_get):
        self.cache = cache if cache is not None else Cache()
        self.http = http

    def text(self, pdb_id: str) -> str:
        pdb_id = validate_pdb_id(pdb_id)
        raw = self.cache.get("rcsb", f"{pdb_id}.pdb")
        if raw is None:
            try:
                raw = self.http(FILE_URL.format(pdb_id=pdb_id))
            except urllib.error.HTTPError as error:
                error.close()
                if error.code == 404:
                    raise ServiceError("structure_not_found", accession=pdb_id) from None
                raise ServiceError("service_unavailable", service=SERVICE) from None
            except (urllib.error.URLError, OSError):
                raise ServiceError("service_unavailable", service=SERVICE) from None
            if b"ATOM" not in raw[:200000] and b"HETATM" not in raw[:200000]:
                raise ServiceError("structure_not_found", accession=pdb_id)
            self.cache.put("rcsb", f"{pdb_id}.pdb", raw)
        return raw.decode("utf-8", errors="replace")

    def structure(self, pdb_id: str) -> Structure:
        return parse_pdb(self.text(pdb_id))
