"""BLAST del EBI: buscar proteínas parecidas en UniProt.

BLAST (Basic Local Alignment Search Tool, 1990) compara una proteína con millones en pocos
segundos: en vez de alinear con todas (como alignment.py, que sería demasiado lento), busca
"semillas" de pocas letras muy parecidas y solo alinea alrededor de ellas.

El EBI lo ofrece como servicio gratuito con trabajos asíncronos (Job Dispatcher):

    1. POST /run                -> identificador del trabajo
    2. GET  /status/<id>        -> QUEUED, RUNNING, FINISHED, ERROR, FAILURE o NOT_FOUND
    3. GET  /result/<id>/json   -> los parecidos encontrados ("hits")

El EBI exige un correo electrónico real en cada búsqueda (para avisar si un uso causa
problemas). Se usa el de Relative Developers (DEFAULT_EMAIL); la variable de entorno
DNAINVENT_EMAIL lo sustituye si hace falta (blast_email()).

Bases de datos:
    uniprotkb_swissprot   ~570.000 proteínas revisadas a mano: rápida y fiable
    uniprotkb             más de 200 millones (incluye las anotadas automáticamente): lenta
"""

import json
import re
import time
import urllib.error
from dataclasses import dataclass
from typing import Callable

from ...errors import InvalidParameterError, ServiceError
from ..cache import Cache, key_for
from ..classify import Hit
from . import http_get, http_post

SERVICE = "BLAST (EBI)"
BASE_URL = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"
SWISSPROT = "uniprotkb_swissprot"
UNIPROTKB = "uniprotkb"
DATABASES = (SWISSPROT, UNIPROTKB)

FIRST_WAIT = 3.0             # s antes de la primera consulta (el EBI pide no insistir)
POLL_EVERY = 5.0             # s entre consultas
TIMEOUT = 15 * 60            # s: UniProtKB completa puede tardar varios minutos
RESULT_MAX_AGE = 30 * 24 * 3600     # UniProt se actualiza cada ~8 semanas

PREFERRED_ORGANISM = "Homo sapiens"
DEFAULT_EMAIL = "relativedevelopers@gmail.com"
EMAIL = re.compile(r"[^@\s]+@[^@\s]+\.[^@\s]+")
# "Hemoglobin subunit beta OS=Homo sapiens OX=9606 GN=HBB PE=1 SV=2"
_UNIPROT_FIELDS = re.compile(r"\s(OS|OX|GN|PE|SV)=")


@dataclass(frozen=True)
class BlastHit:
    """Una proteína real parecida a la del usuario (su mejor tramo alineado)."""

    accession: str           # "P68871"
    entry_name: str          # "HBB_HUMAN"
    description: str         # "Hemoglobin subunit beta"
    organism: str            # "Homo sapiens"
    gene: str                # "HBB"
    length: int              # aminoácidos de la proteína real
    identity: float          # % de aminoácidos idénticos en el tramo alineado
    query_coverage: float    # % de la proteína del usuario dentro del tramo
    subject_coverage: float  # % de la proteína real dentro del tramo
    evalue: float
    bit_score: float
    query_from: int          # tramo alineado (desde 1, como en BLAST)
    query_to: int
    hit_from: int
    hit_to: int
    query_aligned: str = ""  # el tramo con huecos '-', para dibujarlo
    hit_aligned: str = ""
    markup: str = ""         # línea intermedia de BLAST: letra = igual, '+' = parecido

    def to_hit(self) -> Hit:
        return Hit(self.accession, self.identity, self.query_coverage,
                   self.subject_coverage, self.evalue,
                   preferred=self.organism == PREFERRED_ORGANISM)


@dataclass(frozen=True)
class BlastResult:
    database: str
    query: str
    hits: list[BlastHit]


def blast_email() -> str:
    """Correo para el EBI: DNAINVENT_EMAIL si está definida; si no, el de Relative Developers."""
    import os
    return os.environ.get("DNAINVENT_EMAIL") or DEFAULT_EMAIL


def validate_email(email: str | None) -> str:
    email = (email or "").strip()
    if not EMAIL.fullmatch(email):
        raise InvalidParameterError("invalid_email")
    return email


def parse_result(raw: bytes | str, query: str, database: str) -> BlastResult:
    """JSON de resultados del EBI -> BlastResult (solo el mejor tramo de cada proteína)."""
    try:
        data = json.loads(raw)
        hits = [_parse_hit(item, len(query)) for item in data.get("hits", [])]
    except (ValueError, TypeError, KeyError, AttributeError):
        raise ServiceError("search_failed", service=SERVICE) from None
    return BlastResult(database, query, [hit for hit in hits if hit is not None])


def _parse_hit(item: dict, query_length: int) -> BlastHit | None:
    hsps = item.get("hit_hsps") or []
    if not hsps:
        return None
    hsp = hsps[0]                              # el EBI los ordena del mejor al peor
    length = int(item.get("hit_len") or 0)
    q_from, q_to = int(hsp["hsp_query_from"]), int(hsp["hsp_query_to"])
    h_from, h_to = int(hsp["hsp_hit_from"]), int(hsp["hsp_hit_to"])
    fields = _uniprot_fields(item.get("hit_desc") or "")
    return BlastHit(
        accession=str(item.get("hit_acc") or ""),
        entry_name=str(item.get("hit_id") or ""),
        description=item.get("hit_uni_de") or fields.get("DE", ""),
        organism=item.get("hit_uni_os") or fields.get("OS", ""),
        gene=item.get("hit_uni_gn") or fields.get("GN", ""),
        length=length,
        identity=min(100.0, float(hsp["hsp_identity"])),
        query_coverage=_percent(q_to - q_from + 1, query_length),
        subject_coverage=_percent(h_to - h_from + 1, length),
        evalue=float(hsp["hsp_expect"]),
        bit_score=float(hsp.get("hsp_bit_score") or 0.0),
        query_from=q_from, query_to=q_to, hit_from=h_from, hit_to=h_to,
        query_aligned=hsp.get("hsp_qseq") or "",
        hit_aligned=hsp.get("hsp_hseq") or "",
        markup=hsp.get("hsp_mseq") or "",
    )


def _percent(part: int, total: int) -> float:
    return min(100.0, part / total * 100) if total else 0.0


def _uniprot_fields(text: str) -> dict:
    """"Nombre OS=Especie OX=9606 GN=GEN ..." -> {"DE": "Nombre", "OS": "Especie", ...}"""
    parts = _UNIPROT_FIELDS.split(" " + text)
    out = {"DE": parts[0].strip()}
    for key, value in zip(parts[1::2], parts[2::2]):
        out[key] = value.strip()
    return out


class BlastClient:
    def __init__(self, email: str | None, cache: Cache | None = None,
                 get: Callable[[str], bytes] = http_get,
                 post: Callable[[str, dict], bytes] = http_post,
                 sleep: Callable[[float], None] = time.sleep,
                 clock: Callable[[], float] = time.monotonic,
                 timeout: float = TIMEOUT):
        self.email = email
        self.cache = cache if cache is not None else Cache()
        self.get, self.post, self.sleep, self.clock = get, post, sleep, clock
        self.timeout = timeout

    def search(self, protein: str, database: str = SWISSPROT,
               progress: Callable[[str], None] | None = None) -> BlastResult:
        """Busca la proteína en la base de datos (o la lee de la caché).

        progress recibe el estado del trabajo ("QUEUED", "RUNNING") en cada consulta.
        """
        if database not in DATABASES:
            raise ValueError(f"Base de datos desconocida: {database}")
        key = key_for(f"{database}:{protein}") + ".json"
        raw = self.cache.get("blast", key, max_age=RESULT_MAX_AGE)
        if raw is None:
            raw = self._run(protein, database, progress)
            result = parse_result(raw, protein, database)
            self.cache.put("blast", key, raw)
            return result
        return parse_result(raw, protein, database)

    def _run(self, protein: str, database: str, progress) -> bytes:
        email = validate_email(self.email)
        job = self._call(self.post, f"{BASE_URL}/run", {
            "email": email, "program": "blastp", "stype": "protein",
            "sequence": protein, "database": database,
        }).decode().strip()
        started = self.clock()
        self.sleep(FIRST_WAIT)
        while True:
            status = self._call(self.get, f"{BASE_URL}/status/{job}").decode().strip().upper()
            if status == "FINISHED":
                return self._call(self.get, f"{BASE_URL}/result/{job}/json")
            if status not in ("QUEUED", "RUNNING", "PENDING"):
                raise ServiceError("search_failed", service=SERVICE)    # ERROR, FAILURE...
            if progress:
                progress(status)
            if self.clock() - started > self.timeout:
                raise ServiceError("search_timeout", service=SERVICE,
                                   minutes=round(self.timeout / 60))
            self.sleep(POLL_EVERY)

    @staticmethod
    def _call(function, *args) -> bytes:
        try:
            return function(*args)
        except urllib.error.HTTPError as error:
            if 400 <= error.code < 500:        # parámetros rechazados (p. ej. el correo)
                raise ServiceError("search_failed", service=SERVICE) from None
            raise ServiceError("service_unavailable", service=SERVICE) from None
        except (urllib.error.URLError, OSError):
            raise ServiceError("service_unavailable", service=SERVICE) from None
