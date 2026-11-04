"""Cliente común para las E-utilities del NCBI (Gene, Nucleotide, ClinVar).

Normas del NCBI que cumple:
  * Como máximo 3 peticiones por segundo sin API key y 10 con ella. La key es gratuita
    (cuenta del NCBI -> Settings -> API Key Management) y se lee de la variable de entorno
    NCBI_API_KEY: nunca va en el código.
  * Cada petición lleva `tool` y `email` (el mismo correo que BLAST: DNAINVENT_EMAIL o el de
    Relative Developers).
  * Si el NCBI responde 429 (demasiadas peticiones) o 5xx (ocupado), se espera y se reintenta.

Recibe la función de descarga como parámetro (como los clientes de dnainvent.protein), así se
prueba sin red.
"""

import json
import os
import time
import urllib.error
import urllib.parse
from typing import Callable

from ..errors import ServiceError
from ..protein.clients import http_get
from ..protein.clients.blast import blast_email

EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOOL_NAME = "DNAinvent"
MAX_RETRIES = 3
INTERVAL_WITHOUT_KEY = 0.34      # s entre peticiones: 3 por segundo, con margen
INTERVAL_WITH_KEY = 0.11         # 10 por segundo


class EUtilsClient:
    service = "NCBI"             # nombre que ven los mensajes de error

    def __init__(self, http: Callable[[str], bytes] = http_get, email: str | None = None,
                 api_key: str | None = None, sleep: Callable[[float], None] = time.sleep):
        self.http = http
        self.email = email or blast_email()
        self.api_key = os.environ.get("NCBI_API_KEY") if api_key is None else api_key
        self.sleep = sleep
        self.min_interval = INTERVAL_WITH_KEY if self.api_key else INTERVAL_WITHOUT_KEY
        self._last_request = 0.0

    def url(self, endpoint: str, params: dict) -> str:
        params = dict(params, tool=TOOL_NAME, email=self.email)
        if self.api_key:
            params["api_key"] = self.api_key
        return EUTILS_URL + endpoint + "?" + urllib.parse.urlencode(params)

    def request(self, endpoint: str, params: dict) -> bytes:
        url = self.url(endpoint, params)
        for attempt in range(1, MAX_RETRIES + 1):
            wait = self.min_interval - (time.monotonic() - self._last_request)
            if wait > 0:
                self.sleep(wait)
            self._last_request = time.monotonic()
            try:
                return self.http(url)
            except urllib.error.HTTPError as error:
                error.close()
                retry = error.code == 429 or error.code >= 500
            except (urllib.error.URLError, OSError):      # sin red, DNS, tiempo agotado...
                retry = True
            if not retry or attempt == MAX_RETRIES:
                break
            self.sleep(2 ** attempt)
        raise ServiceError("service_unavailable", service=self.service)

    def json(self, endpoint: str, params: dict) -> dict:
        try:
            data = json.loads(self.request(endpoint, dict(params, retmode="json")))
        except ValueError:
            raise ServiceError("service_unavailable", service=self.service) from None
        if not isinstance(data, dict) or "error" in data:
            raise ServiceError("search_failed", service=self.service)
        return data
