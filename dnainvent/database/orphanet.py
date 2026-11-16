"""Orphanet (Orphadata): cómo se hereda una enfermedad rara y lo frecuente que es.

    provider = OrphanetProvider()
    stats = provider.apply(conn)          # rellena diseases.inheritance, onset y prevalence_*

Datos: Orphanet, «Orphadata», licencia CC BY 4.0 (https://www.orphadata.com/). Hay que citarlo
donde se muestren. Solo cubre ENFERMEDADES RARAS: el Alzheimer o la intolerancia a la lactosa,
por ejemplo, no están, y se quedan sin estos datos.

No es una API: son dos archivos XML que se descargan enteros (unos pocos MB) y se guardan en la
caché 30 días (~/.dnainvent/cache/orphanet):

    en_product9_ages.xml   historia natural: tipo de herencia y edad de inicio
    en_product9_prev.xml   epidemiología: prevalencia, por tipo, zona y fuente

Estructura (descrita en «Orphadata Science Data Description»):

    <JDBOR date="2026-06-23 07:57:58" ...>
      <DisorderList count="...">
        <Disorder id="...">
          <OrphaCode>586</OrphaCode>
          <Name lang="en">Cystic fibrosis</Name>
          ... <TypeOfInheritance><Name lang="en">Autosomal recessive</Name></TypeOfInheritance>
          ... <AverageAgeOfOnset><Name lang="en">Infancy</Name></AverageAgeOfOnset>
          <PrevalenceList count="2">
            <Prevalence id="...">
              <Source>...</Source>
              <PrevalenceType><Name lang="en">Point prevalence</Name></PrevalenceType>
              <PrevalenceQualification><Name lang="en">Value and class</Name></...>
              <PrevalenceClass><Name lang="en">1-5 / 10 000</Name></PrevalenceClass>
              <ValMoy>1.2</ValMoy>
              <PrevalenceGeographic><Name lang="en">Europe</Name></PrevalenceGeographic>
              <PrevalenceValidationStatus><Name lang="en">Validated</Name></...>

Se lee buscando las etiquetas por su nombre, estén donde estén dentro de cada Disorder: así no
se rompe si Orphanet añade un nivel intermedio (TypeOfInheritanceList...).

Las enfermedades se enlazan por el código Orphanet que da ClinVar (diseases.orpha_codes). Si
una enfermedad no lo tiene, se prueba con el nombre exacto (sin mayúsculas).
"""

import sqlite3
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

from ..errors import ServiceError
from ..protein.cache import Cache
from ..protein.clients import http_get
from . import get_or_create_source

SERVICE = "Orphanet"
MAX_AGE = 30 * 24 * 3600
TIMEOUT = 120                 # s: los archivos ocupan varios MB
FILES = {"ages": "en_product9_ages.xml", "prev": "en_product9_prev.xml"}
GITHUB_FOLDERS = {"ages": "Rare%20diseases%20natural%20history",
                  "prev": "Rare%20disease%20epidemiology"}
# Primero la web de Orphadata; si falla, la copia oficial en GitHub (Orphanet/Orphadata_aggregated,
# que guarda los XML con Git LFS: por eso media.githubusercontent.com y no raw.githubusercontent.com)
URLS = {kind: (f"https://www.orphadata.com/data/xml/{name}",
               "https://media.githubusercontent.com/media/Orphanet/Orphadata_aggregated/master/"
               f"Epidemiological%20data/{GITHUB_FOLDERS[kind]}/{name}")
        for kind, name in FILES.items()}

# Valores que en realidad dicen «no se sabe»
NO_DATA = {"", "no data available", "not yet documented", "not applicable", "unknown",
           "unknown_epidemiological_range", "no age of onset data available"}

PREVALENCE_TYPES = ("point prevalence", "birth prevalence", "annual incidence",
                    "lifetime prevalence", "lifelong prevalence")
PLACES = ("worldwide", "europe")

Progress = Callable[[str, dict], None]


def _no_progress(key: str, data: dict) -> None:
    pass


def _slow_http_get(url: str) -> bytes:
    return http_get(url, timeout=TIMEOUT)


def _text(element, path: str) -> str:
    return (element.findtext(path) or "").strip()


def _known(value: str) -> bool:
    return value.lower() not in NO_DATA


def _parse(raw: bytes) -> ET.Element:
    try:
        return ET.fromstring(raw)
    except ET.ParseError:
        raise ServiceError("service_unavailable", service=SERVICE) from None


def _disorders(root: ET.Element):
    """(código Orphanet, nombre, elemento) de cada enfermedad del archivo."""
    for disorder in root.iter("Disorder"):
        code = _text(disorder, "OrphaCode")
        if code:
            yield code, _text(disorder, "Name"), disorder


@dataclass(frozen=True)
class Prevalence:
    type: str                 # 'Point prevalence'
    qualification: str        # 'Value and class', 'Class only', 'Case', 'Family'
    prevalence_class: str     # '1-5 / 10 000'
    value: str                # ValMoy: media, tal cual (puede ser '0.0')
    geographic: str           # 'Worldwide', 'Europe', 'Spain'...
    status: str               # 'Validated', 'Not yet validated'
    source: str


@dataclass
class NaturalHistory:
    name: str = ""
    inheritance: list[str] = field(default_factory=list)
    onset: list[str] = field(default_factory=list)


def parse_natural_history(raw: bytes) -> tuple[dict[str, NaturalHistory], str | None]:
    """en_product9_ages.xml -> ({código: herencia y edad de inicio}, fecha del archivo)."""
    root = _parse(raw)
    out: dict[str, NaturalHistory] = {}
    for code, name, disorder in _disorders(root):
        entry = NaturalHistory(name=name)
        for tag, target in (("TypeOfInheritance", entry.inheritance),
                            ("AverageAgeOfOnset", entry.onset)):
            for item in disorder.iter(tag):
                value = _text(item, "Name")
                if _known(value) and value not in target:
                    target.append(value)
        out[code] = entry
    return out, root.get("date")


def parse_prevalence(raw: bytes) -> tuple[dict[str, tuple[str, list[Prevalence]]], str | None]:
    """en_product9_prev.xml -> ({código: (nombre, [prevalencias])}, fecha del archivo)."""
    root = _parse(raw)
    out: dict[str, tuple[str, list[Prevalence]]] = {}
    for code, name, disorder in _disorders(root):
        items = [Prevalence(
            type=_text(p, "PrevalenceType/Name"),
            qualification=_text(p, "PrevalenceQualification/Name"),
            prevalence_class=_text(p, "PrevalenceClass/Name"),
            value=_text(p, "ValMoy"),
            geographic=_text(p, "PrevalenceGeographic/Name"),
            status=_text(p, "PrevalenceValidationStatus/Name"),
            source=_text(p, "Source"),
        ) for p in disorder.iter("Prevalence")]
        out[code] = (name, items)
    return out, root.get("date")


def best_prevalence(items: list[Prevalence]) -> Prevalence | None:
    """La prevalencia que se enseña: con una clase conocida (1-5 / 10 000...), validada, mundial
    o europea antes que de un país, y de tipo «prevalencia puntual» antes que de nacimiento o
    incidencia. Las que solo cuentan casos o familias no dan una proporción y no sirven."""
    usable = [p for p in items if _known(p.prevalence_class)]
    if not usable:
        return None

    def rank(p: Prevalence):
        kind = p.type.lower()
        place = p.geographic.lower()
        return (p.status.lower() != "validated",
                PLACES.index(place) if place in PLACES else len(PLACES),
                PREVALENCE_TYPES.index(kind) if kind in PREVALENCE_TYPES else len(PREVALENCE_TYPES))

    return min(usable, key=rank)


@dataclass
class OrphanetStats:
    diseases: int = 0          # enfermedades de la base de datos
    by_code: int = 0           # enlazadas por código Orphanet
    by_name: int = 0           # enlazadas por nombre
    inheritance: int = 0       # con tipo de herencia
    prevalence: int = 0        # con prevalencia
    file_dates: dict = field(default_factory=dict)   # {'ages': '2026-06-23 ...', 'prev': ...}


class OrphanetProvider:
    service = SERVICE

    def __init__(self, http: Callable[[str], bytes] = _slow_http_get, cache: Cache | None = None,
                 local_dir=None):
        self.http = http
        self.cache = cache if cache is not None else Cache()
        self.local_dir = local_dir      # carpeta con los XML ya descargados a mano (opcional)

    def fetch(self, kind: str, progress: Progress = _no_progress) -> bytes:
        name = FILES[kind]
        if self.local_dir is not None:
            path = Path(self.local_dir) / name
            if path.is_file():
                return path.read_bytes()
        raw = self.cache.get("orphanet", name, max_age=MAX_AGE)
        if raw is not None:
            return raw
        progress("database.stage.orphanet", {"file": name})
        for url in URLS[kind]:
            try:
                raw = self.http(url)
            except OSError:              # HTTPError, URLError y tiempo agotado son OSError
                continue
            if b"<JDBOR" in raw[:5000]:      # y no una página de error o un puntero de Git LFS
                self.cache.put("orphanet", name, raw)
                return raw
        raise ServiceError("service_unavailable", service=SERVICE)

    def apply(self, conn: sqlite3.Connection, progress: Progress = _no_progress) -> OrphanetStats:
        """Rellena la herencia, la edad de inicio y la prevalencia de las enfermedades guardadas."""
        history, ages_date = parse_natural_history(self.fetch("ages", progress))
        prevalence, prev_date = parse_prevalence(self.fetch("prev", progress))
        stats = OrphanetStats(file_dates={"ages": ages_date, "prev": prev_date})
        by_name: dict[str, str] = {}
        for code, entry in history.items():
            by_name.setdefault(entry.name.lower(), code)
        for code, (name, _) in prevalence.items():
            by_name.setdefault(name.lower(), code)

        with conn:
            get_or_create_source(conn, "orphanet")
            for disease in conn.execute("SELECT id, name_en, orpha_codes FROM diseases").fetchall():
                stats.diseases += 1
                codes = [c for c in (disease["orpha_codes"] or "").split("|") if c]
                if codes:
                    stats.by_code += any(c in history or c in prevalence for c in codes)
                elif disease["name_en"] and disease["name_en"].lower() in by_name:
                    codes = [by_name[disease["name_en"].lower()]]
                    stats.by_name += 1
                    conn.execute("UPDATE diseases SET orpha_codes = ? WHERE id = ?",
                                 (codes[0], disease["id"]))
                if not codes:
                    continue
                # El primer código que tenga cada dato (ClinVar pone primero el principal)
                entry = next((history[c] for c in codes if c in history and
                              (history[c].inheritance or history[c].onset)), None)
                best = next((b for b in (best_prevalence(prevalence[c][1]) for c in codes
                                         if c in prevalence) if b), None)
                conn.execute(
                    """UPDATE diseases SET inheritance = ?, onset = ?, prevalence_class = ?,
                       prevalence_type = ?, prevalence_geo = ? WHERE id = ?""",
                    ("|".join(entry.inheritance) or None if entry else None,
                     "|".join(entry.onset) or None if entry else None,
                     best.prevalence_class if best else None, best.type if best else None,
                     best.geographic if best else None, disease["id"]))
                stats.inheritance += bool(entry and entry.inheritance)
                stats.prevalence += best is not None
        return stats
