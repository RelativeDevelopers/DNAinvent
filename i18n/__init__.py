"""Traducciones de DNAinvent.

Cada idioma es un archivo JSON en i18n/locales/ con las mismas claves.
Las claves son anidadas y se piden con puntos:

    tr = Translator("en")
    tr("translate.start_at", position=4)   ->  "Start at position 4."

Los parámetros se escriben entre llaves en el JSON: "Start at position {position}."
Es una sintaxis sencilla a propósito, para que la web (JavaScript) pueda usar
exactamente los mismos archivos.

Idiomas: español (es, el de referencia), inglés (en) y catalán (ca). Para añadir otro basta
con crear su JSON: los tests comprueban que tenga todas las claves.
"""

import json
import re
from functools import lru_cache
from pathlib import Path

from dnainvent.errors import DNAinventError

LOCALES_DIR = Path(__file__).parent / "locales"
DEFAULT_LANGUAGE = "es"
PLACEHOLDER = re.compile(r"\{(\w+)\}")
_SUPERSCRIPTS = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def available_languages() -> dict[str, str]:
    """Código -> nombre del idioma, p. ej. {"en": "English", "es": "Español"}."""
    return {path.stem: _load(path.stem)["_meta"]["name"]
            for path in sorted(LOCALES_DIR.glob("*.json"))}


@lru_cache
def _load(language: str) -> dict:
    path = LOCALES_DIR / f"{language}.json"
    if not path.exists():
        raise ValueError(f"Idioma no disponible: {language}")
    return json.loads(path.read_text(encoding="utf-8"))


def _lookup(data: dict, key: str):
    node = data
    for part in key.split("."):
        if not isinstance(node, dict) or part not in node:
            return None
        node = node[part]
    return node


class Translator:
    def __init__(self, language: str = DEFAULT_LANGUAGE):
        self._fallback = _load(DEFAULT_LANGUAGE)
        self.set_language(language)

    def set_language(self, language: str) -> None:
        self._data = _load(language)        # ValueError si no existe
        self.language = language

    def next_language(self) -> str:
        """Pasa al siguiente idioma (tecla L del microscopio) y lo devuelve."""
        codes = list(available_languages())
        self.set_language(codes[(codes.index(self.language) + 1) % len(codes)])
        return self.language

    def __call__(self, key: str, **params) -> str:
        """Texto de la clave con los parámetros sustituidos.
        Si falta en este idioma usa el idioma por defecto, y si tampoco existe, la propia clave."""
        text = _lookup(self._data, key)
        if not isinstance(text, str):
            text = _lookup(self._fallback, key)
        if not isinstance(text, str):
            return key
        return self._fill(text, params)

    def lines(self, key: str, **params) -> list[str]:
        """Para las claves que son una lista de líneas (las explicaciones del microscopio)."""
        value = _lookup(self._data, key)
        if not isinstance(value, list):
            value = _lookup(self._fallback, key)
        if not isinstance(value, list):
            return [key]
        return [self._fill(line, params) for line in value]

    def _fill(self, text: str, params: dict) -> str:
        return PLACEHOLDER.sub(lambda m: self._format(params.get(m.group(1), m.group(0))), text)

    def _format(self, value) -> str:
        if isinstance(value, bool):
            return self("common.yes" if value else "common.no")
        if isinstance(value, int):
            return self.number(value)
        if isinstance(value, (list, tuple)):
            return ", ".join(str(item) for item in value)
        return str(value)

    def number(self, value: int) -> str:
        """Entero con separador de miles del idioma: 12.345 (es) o 12,345 (en)."""
        return f"{value:,}".replace(",", self._data["_meta"]["thousands_separator"])

    def decimal(self, value: float, digits: int = 2, trim: bool = False) -> str:
        """Decimal con el separador del idioma: 3,25 (es) o 3.25 (en).
        Con trim=True quita los ceros sobrantes: 1,750 -> 1,75 y 2,000 -> 2."""
        text = f"{value:.{digits}f}"
        if trim and "." in text:
            text = text.rstrip("0").rstrip(".")
        return text.replace(".", self._data["_meta"]["decimal_separator"])

    @property
    def decimal_separator(self) -> str:
        return self._data["_meta"]["decimal_separator"]

    def sci(self, value: float, digits: int = 1) -> str:
        """Notación científica legible: 1.2e-10 -> '1,2×10⁻¹⁰' (es) o '1.2×10⁻¹⁰' (en)."""
        if value == 0:
            return "0"
        mantissa, exponent = f"{value:.{digits}e}".split("e")
        return (mantissa.replace(".", self.decimal_separator) + "×10"
                + str(int(exponent)).translate(_SUPERSCRIPTS))

    def count(self, n: int, noun: str) -> str:
        """'1 protón', '6 protones': las dos formas están en units.<noun>.one / .other."""
        return self(f"units.{noun}.{'one' if n == 1 else 'other'}", n=n)

    # --- nombres de los elementos (vienen de AtomicWorld, en en/es/ca) ---------------
    def element(self, number: int) -> str:
        """'Carbono', 'Carbon', 'Carboni'. Solo lo usa el microscopio (necesita atomicworld)."""
        from atomicworld import names
        language = self.language if self.language in names.LANGUAGES else "en"
        return names.element(number, language)

    def of(self, name: str) -> str:
        """'de carbono' (es), "d'oxigen" / 'de carboni' (ca), 'Carbon' (en, sin preposición).

        En castellano y catalán los nombres de los elementos van en minúscula dentro de la frase.
        """
        if self.language == "es":
            return f"de {name.lower()}"
        if self.language == "ca":
            name = name.lower()
            return f"d'{name}" if _catalan_elides(name) else f"de {name}"
        return name

    def error(self, error: DNAinventError) -> str:
        """Mensaje de un error del núcleo en este idioma."""
        return self.error_message(error.code, error.params)

    def error_message(self, code: str, params: dict) -> str:
        """Mensaje de cualquier error a partir de su código y sus datos.
        Es el mismo formato que devuelve la API: {"code": ..., "params": {...}}.
        'kind' y 'field' son nombres internos, así que también se traducen."""
        params = dict(params)
        if "kind" in params:
            params["kind"] = self(f"kinds.{params['kind']}")
        if "field" in params:
            params["field"] = self(f"fields.{params['field']}")
        return self(f"errors.{code}", **params)

    def amino_acid(self, letter: str) -> str:
        """Nombre del aminoácido con su código de 3 letras: "Glutamina (Gln)"."""
        from dnainvent.constants import AMINO_ACID_CODES
        return f"{self(f'amino_acids.{letter}')} ({AMINO_ACID_CODES[letter]})"


_VOWELS = "aeiouàèéíïòóúü"


def _catalan_elides(word: str) -> bool:
    """En catalán 'de' se apostrofa ante vocal o h + vocal (d'oxigen, d'hidrogen), salvo ante
    i/u seguidas de vocal, que suenan como consonante (de iode)."""
    w = word.lower()
    if w.startswith("h"):
        w = w[1:]
    if not w or w[0] not in _VOWELS:
        return False
    return not (w[0] in "iu" and len(w) > 1 and w[1] in _VOWELS)
