"""ESMFold (Meta, 2022): predice la estructura a partir de una sola secuencia.

A diferencia de AlphaFold, no necesita buscar proteínas parecidas (alineamiento múltiple):
un modelo de lenguaje de proteínas (ESM-2) ha "aprendido" la evolución leyendo millones de
secuencias. Por eso es rápido (segundos) y sirve para proteínas inventadas, aunque es algo
menos preciso que AlphaFold.

    POST https://api.esmatlas.com/foldSequence/v1/pdb/      cuerpo: la secuencia, tal cual
    -> archivo PDB con el pLDDT en la columna B

El servidor público limita la longitud y el número de peticiones, y alguna vez ha estado
caído: los errores se convierten en ServiceError y el resto del programa sigue.

Lin et al. 2023, Science 379: 1123-1130.
"""

import urllib.error
from typing import Callable

from ...errors import InvalidParameterError, ServiceError
from ...validation import validate_protein
from ..cache import Cache, key_for
from ..clients import http_post_text
from ..pdb import Atom, Structure, parse_pdb

API_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
MAX_LENGTH = 400          # aminoácidos: el límite que documentan los ejemplos del servidor
TIMEOUT = 120             # s: una proteína de 400 aminoácidos tarda unos segundos


class ESMFold:
    name = "ESMFold"
    max_length = MAX_LENGTH

    def __init__(self, cache: Cache | None = None,
                 post: Callable[[str, str], bytes] | None = None):
        self.cache = cache if cache is not None else Cache()
        self.post = post or (lambda url, text: http_post_text(url, text, timeout=TIMEOUT))

    def predict(self, protein: str) -> Structure:
        protein = validate_protein(protein)
        if len(protein) > self.max_length:
            raise InvalidParameterError("prediction_too_long", length=len(protein),
                                        maximum=self.max_length)
        key = key_for(protein) + ".pdb"
        raw = self.cache.get("esmfold", key)
        if raw is None:
            raw = self._download(protein)
        structure = _normalize(parse_pdb(raw.decode("utf-8", errors="replace")))
        if structure.sequence() != protein:          # respuesta vacía, de error o truncada
            raise ServiceError("prediction_failed", service=self.name)
        self.cache.put("esmfold", key, raw)
        return structure

    def _download(self, protein: str) -> bytes:
        try:
            return self.post(API_URL, protein)
        except urllib.error.HTTPError as error:
            # 4xx: el servidor responde pero rechaza esta proteína; 5xx: el servidor falla
            code = "prediction_failed" if 400 <= error.code < 500 else "service_unavailable"
            raise ServiceError(code, service=self.name) from None
        except (urllib.error.URLError, OSError):     # sin red, SSL, tiempo agotado...
            raise ServiceError("service_unavailable", service=self.name) from None


def _normalize(structure: Structure) -> Structure:
    """pLDDT siempre de 0 a 100. Algunas versiones de ESMFold lo escriben de 0 a 1."""
    if not structure.atoms or max(a.bfactor for a in structure.atoms) > 1.0:
        title = structure.title or "ESMFOLD PREDICTION"
        return Structure(structure.atoms, title)
    atoms = [_scaled(a) for a in structure.atoms]
    return Structure(atoms, structure.title or "ESMFOLD PREDICTION")


def _scaled(atom: Atom) -> Atom:
    from dataclasses import replace
    return replace(atom, bfactor=round(atom.bfactor * 100, 2))
