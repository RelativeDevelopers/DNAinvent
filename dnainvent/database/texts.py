"""Textos escritos a mano para los genes de serie: nombres en cada idioma, explicaciones,
«¿Sabías que…?» y «por qué funciona» de los fármacos.

    from dnainvent.database import texts
    texts.gene("HBB", "es")["name"]           # 'Hemoglobina, subunidad beta'
    texts.disease_name("Hb SS disease", "ca") # 'Anèmia falciforme (drepanocitosi)'
    texts.drug("CHEMBL4101807", "en")["why"]  # por qué funciona
    texts.sources("drugs", "CHEMBL4101807")   # ['https://...']
    texts.disease_note("G6PD deficiency", "es")  # qué se hace cuando no hay fármacos

Están en dnainvent/database/texts/:
    es.json, ca.json, en.json   los textos de cada idioma (mismas claves en los tres)
    common.json                 las fuentes de cada texto (iguales en todos los idiomas)

Claves:
    genes       símbolo ('HBB')
    variants    'SÍMBOLO:c.HGVS' ('HBB:c.20A>T')
    diseases    nombre en inglés tal como lo da ClinVar ('Hb SS disease')
    drugs       ID de ChEMBL ('CHEMBL4101807')
    disease_notes  nombre en inglés de ClinVar: una nota sobre el tratamiento de la enfermedad
    shared      textos de «por qué funciona» que comparten varios fármacos (los concentrados
                de factor VIII): en common.json el fármaco lleva {"same_as": "factor_viii"}

Si un texto no está, se devuelve None y las fichas siguen con lo que dan las fuentes
(en inglés, con la marca «(nombre en inglés)»). Los JSON se pueden corregir sin tocar código.
"""

import json
from functools import lru_cache
from pathlib import Path

TEXTS_DIR = Path(__file__).parent / "texts"
LANGUAGES = ("es", "ca", "en")


@lru_cache(maxsize=None)
def _load(name: str) -> dict:
    path = TEXTS_DIR / f"{name}.json"
    if not path.is_file():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _get(section: str, key: str, language: str):
    return _load(language).get(section, {}).get(key) if language in LANGUAGES else None


def gene(symbol: str, language: str) -> dict | None:
    """{'name': ..., 'summary': ..., 'did_you_know': [...]}"""
    return _get("genes", symbol.upper(), language)


def variant(symbol: str, hgvs_c: str | None, language: str) -> str | None:
    return _get("variants", f"{symbol.upper()}:{hgvs_c}", language) if hgvs_c else None


def disease_name(name_en: str | None, language: str) -> str | None:
    return _get("diseases", name_en, language) if name_en else None


def _shared(section: str, key: str) -> str | None:
    """La clave compartida de un fármaco (common.json: {"same_as": ...}), si la tiene."""
    return (_load("common").get(section, {}).get(key) or {}).get("same_as")


def drug(chembl_id: str | None, language: str) -> dict | None:
    """{'name': ..., 'why': ...}"""
    entry = _get("drugs", chembl_id, language) if chembl_id else None
    shared = _shared("drugs", chembl_id) if entry else None
    if shared and not entry.get("why"):
        entry = {**entry, "why": _get("shared", shared, language)}
    return entry


def disease_note(name_en: str | None, language: str) -> str | None:
    return _get("disease_notes", name_en, language) if name_en else None


def sources(section: str, key: str) -> list[str]:
    common = _load("common")
    entry = common.get(section, {}).get(key) or {}
    if entry.get("same_as"):
        entry = common.get("shared", {}).get(entry["same_as"]) or {}
    return list(entry.get("sources", []))
