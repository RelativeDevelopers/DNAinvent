"""Los nombres escritos a mano de las partes del cuerpo y de los tipos de célula.

    from dnainvent.anatomy import names
    names.concept("FMA7088", "es")        # 'Corazón' (None si aún no está traducido)
    names.cell("CL:0000746", "ca")        # 'cèl·lula del múscul cardíac'
    names.normalize("Corazón")            # 'corazon' (para buscar sin acentos)

Están en dnainvent/anatomy/texts/{es,ca,en}.json, por el ID de la FMA (partes del cuerpo) y de
la Cell Ontology (células). BodyParts3D y las tablas ASCT+B dan los nombres en inglés: lo que no
está aquí se enseña en inglés, con la marca «(nombre en inglés)», como en «Genes reales».
"""

import json
import unicodedata
from functools import lru_cache
from pathlib import Path

TEXTS_DIR = Path(__file__).parent / "texts"


@lru_cache(maxsize=None)
def _load(language: str) -> dict:
    path = TEXTS_DIR / f"{language}.json"
    if not path.is_file():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def concept(fma: str, language: str) -> str | None:
    return _load(language).get("concepts", {}).get(fma)


def cell(cl: str, language: str) -> str | None:
    return _load(language).get("cells", {}).get(cl)


def concepts(language: str) -> dict[str, str]:
    return dict(_load(language).get("concepts", {}))


def english(name: str) -> str:
    """BodyParts3D escribe casi todo en minúsculas: «left femur» -> «Left femur»."""
    return name[:1].upper() + name[1:] if name else name


def normalize(text: str) -> str:
    """Minúsculas y sin acentos, para comparar («Fémur» = «femur»)."""
    text = unicodedata.normalize("NFKD", text or "")
    return "".join(c for c in text if not unicodedata.combining(c)).lower().replace("·", "").strip()
