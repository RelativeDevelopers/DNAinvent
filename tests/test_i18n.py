"""Tests de idiomas: que todos los archivos estén completos y sean coherentes entre sí."""

import json
import re
import tempfile
from pathlib import Path

import pytest

from dnainvent import DNAinventError, alignment, bits, crispr, helix3d, sequences, tools
from dnainvent.constants import AMINO_ACID_CODES
from dnainvent.errors import ERROR_CODES
from dnainvent.protein import orf
from dnainvent.protein.cache import Cache
from dnainvent.protein.clients.afdb import AlphaFoldDB
from dnainvent.protein.predictors.esmfold import ESMFold
from i18n import LOCALES_DIR, Translator, available_languages

LANGUAGES = sorted(available_languages())
PLACEHOLDER = re.compile(r"\{(\w+)\}")


def flatten(data: dict, prefix: str = "") -> dict[str, str]:
    """{"a": {"b": "x"}} -> {"a.b": "x"}"""
    out = {}
    for key, value in data.items():
        full = f"{prefix}{key}"
        if isinstance(value, dict):
            out.update(flatten(value, full + "."))
        else:
            out[full] = value
    return out


def load(language: str) -> dict[str, str]:
    return flatten(json.loads((LOCALES_DIR / f"{language}.json").read_text(encoding="utf-8")))


def test_there_are_at_least_two_languages():
    assert {"es", "en"} <= set(LANGUAGES)


@pytest.mark.parametrize("language", LANGUAGES)
def test_same_keys_in_every_language(language):
    reference = load("es")
    assert set(load(language)) == set(reference)


@pytest.mark.parametrize("language", LANGUAGES)
def test_same_placeholders_in_every_language(language):
    # Si en español una frase usa {count}, en inglés también debe usarlo:
    # si no, el número desaparecería o la frase mostraría "{count}" tal cual.
    reference = load("es")
    for key, text in load(language).items():
        if isinstance(text, str):
            assert set(PLACEHOLDER.findall(text)) == set(PLACEHOLDER.findall(reference[key])), key


@pytest.mark.parametrize("language", LANGUAGES)
def test_every_error_code_is_translated(language):
    keys = load(language)
    for code in ERROR_CODES:
        assert f"errors.{code}" in keys, code


@pytest.mark.parametrize("language", LANGUAGES)
def test_every_amino_acid_is_translated(language):
    keys = load(language)
    for letter in AMINO_ACID_CODES:
        assert f"amino_acids.{letter}" in keys, letter


def test_every_key_used_in_console_exists():
    console = Path(__file__).parents[1] / "ui" / "console.py"
    used = set(re.findall(r'"((?:[a-z_]+\.)+[a-z_]+)"', console.read_text(encoding="utf-8")))
    known = set(load("es"))
    missing = {key for key in used if key not in known and not key.endswith(".py")}
    assert not missing, missing


# Cada error real del núcleo, provocado a propósito, y su mensaje en cada idioma.
ERROR_CASES = [
    lambda: sequences.complement(""),
    lambda: sequences.complement("ATX"),
    lambda: sequences.from_binary("001"),
    lambda: crispr.repair_hdr("ACGT", 10, "A"),
    lambda: crispr.find_targets("ACGT", 2),
    lambda: bits.generate_bits(0),
    lambda: tools.count_permutations("ACG", -1),
    lambda: helix3d.build_helix("ACGT", start=9),
    lambda: orf.best_orf("CCCC"),
    lambda: orf.find_orfs("ATG", min_length=0),
    lambda: alignment.align("MK1", "MK", kind="protein"),
    lambda: AlphaFoldDB(Cache(tempfile.mkdtemp()), http=lambda url: None).entry("HBB"),
    lambda: AlphaFoldDB(Cache(tempfile.mkdtemp()), http=lambda url: b"[]").entry("P00000"),
    lambda: AlphaFoldDB(Cache(tempfile.mkdtemp()), http=lambda url: b"no es json").entry("P00001"),
    lambda: ESMFold(Cache(tempfile.mkdtemp())).predict("M" * 401),
    lambda: ESMFold(Cache(tempfile.mkdtemp()), post=lambda url, text: b"").predict("MKV"),
]


@pytest.mark.parametrize("language", LANGUAGES)
@pytest.mark.parametrize("case", ERROR_CASES)
def test_core_errors_produce_complete_messages(language, case):
    with pytest.raises(DNAinventError) as info:
        case()
    message = Translator(language).error(info.value)
    assert not message.startswith("errors.")          # la clave existe
    assert not PLACEHOLDER.search(message)            # no queda ningún {hueco}


def test_errors_can_be_sent_as_json():
    with pytest.raises(DNAinventError) as info:
        sequences.complement("ATX")
    assert json.loads(json.dumps(info.value.to_dict())) == {
        "code": "invalid_characters",
        "params": {"kind": "dna", "chars": ["X"], "allowed": ["A", "C", "G", "T"]},
    }


def test_unknown_error_code_is_rejected():
    with pytest.raises(ValueError):
        DNAinventError("this_code_does_not_exist")


def test_number_formatting_per_language():
    assert Translator("es").number(1234567) == "1.234.567"
    assert Translator("en").number(1234567) == "1,234,567"
    assert Translator("es").decimal(90.909, 2) == "90,91"
    assert Translator("en").decimal(1.75, 3, trim=True) == "1.75"
    assert Translator("es").decimal(2.0, 3, trim=True) == "2"


def test_missing_key_falls_back_to_key():
    assert Translator("en")("does.not.exist") == "does.not.exist"


def test_amino_acid_name():
    assert Translator("es").amino_acid("Q") == "Glutamina (Gln)"
    assert Translator("en").amino_acid("Q") == "Glutamine (Gln)"


@pytest.mark.parametrize("language", LANGUAGES)
def test_every_confidence_band_is_translated(language):
    from dnainvent.protein import confidence
    keys = load(language)
    for _, band, _ in confidence.BANDS:
        assert f"protein.confidence.{band}" in keys, band
