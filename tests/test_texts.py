"""Tests de los textos escritos a mano (dnainvent/database/texts/*.json)."""

import json
import re

import pytest

from dnainvent.database import texts

LANGS = texts.LANGUAGES


def load(name):
    return json.loads((texts.TEXTS_DIR / f"{name}.json").read_text(encoding="utf-8"))


def keys(data, prefix=""):
    """Todas las claves de un JSON anidado, sin las listas: {'genes.HBB.name', ...}"""
    out = set()
    for key, value in data.items():
        if key.startswith("_"):
            continue
        full = f"{prefix}{key}"
        if isinstance(value, dict):
            out |= keys(value, full + ".")
        else:
            out.add(full)
    return out


def test_same_texts_in_every_language():
    reference = keys(load("es"))
    for lang in LANGS:
        assert keys(load(lang)) == reference, lang


def test_did_you_know_has_the_same_number_of_facts():
    for symbol, gene in load("es")["genes"].items():
        for lang in LANGS:
            assert len(load(lang)["genes"][symbol].get("did_you_know", [])) == len(gene.get("did_you_know", []))


def test_every_text_with_facts_has_sources():
    common = load("common")
    es = load("es")
    for section in ("genes", "drugs", "disease_notes", "shared"):
        for key in es[section]:
            assert texts.sources(section, key), f"{section}.{key} sin fuentes"
    for section, entries in common.items():
        if section.startswith("_"):
            continue
        for key, entry in entries.items():
            if "same_as" in entry:            # usa las fuentes del texto compartido
                assert entry["same_as"] in common["shared"], key
                continue
            assert all(url.startswith("https://") for url in entry["sources"]), key


def test_keys_have_the_right_shape():
    es = load("es")
    assert all(re.fullmatch(r"[A-Z0-9]+", k) for k in es["genes"])
    assert all(re.fullmatch(r"[A-Z0-9]+:c\.\S+", k) for k in es["variants"])
    assert all(re.fullmatch(r"CHEMBL\d+", k) for k in es["drugs"])


def test_no_empty_or_untranslated_placeholders():
    for lang in LANGS:
        text = json.dumps(load(lang), ensure_ascii=False)
        assert '""' not in text and "TODO" not in text
        for section in ("genes", "drugs"):
            for key, entry in load(lang)[section].items():
                assert entry["name"].strip(), (lang, key)


def test_lookups():
    assert texts.gene("hbb", "es")["name"] == "Hemoglobina, subunidad beta"
    assert texts.disease_name("Hb SS disease", "ca") == "Anèmia falciforme (drepanocitosi)"
    assert texts.drug("CHEMBL4101807", "en")["name"] == "Voxelotor"
    assert texts.variant("HBB", "c.20A>T", "es").startswith("La variante de la anemia falciforme")
    assert texts.gene("NOEXISTE", "es") is None and texts.drug(None, "es") is None
    assert texts.gene("HBB", "fr") is None
    # fármacos con un texto compartido (los concentrados de factor VIII)
    octocog = texts.drug("CHEMBL4297873", "es")
    assert octocog["name"] == "Octocog alfa" and octocog["why"].startswith("Es el propio factor VIII")
    assert texts.sources("drugs", "CHEMBL4297873") == texts.sources("shared", "factor_viii")
    assert texts.disease_note("G6PD deficiency", "en").startswith("There are no drugs")


def test_every_gene_of_the_catalog_has_texts():
    from dnainvent.database.catalog import CATALOG_GENES
    for symbol in CATALOG_GENES:
        for lang in LANGS:
            gene = texts.gene(symbol, lang)
            assert gene and gene["summary"] and gene["did_you_know"], (symbol, lang)


def test_console_uses_the_texts(tmp_path, monkeypatch, capsys):
    from test_catalog import make_builder
    from test_database import FakeNCBI
    from dnainvent.database import open_database
    from i18n import Translator
    from ui import console
    conn = open_database(tmp_path / "catalog.db")
    make_builder(conn, FakeNCBI(), tmp_path).run(["HBB"])
    conn.close()
    monkeypatch.setenv("DNAINVENT_CATALOG", str(tmp_path / "catalog.db"))
    monkeypatch.setattr(console, "_catalog", None)
    monkeypatch.setattr(console, "tr", Translator("ca"))
    answers = iter(["1", "1", "HbS", "0", "0", ""])
    monkeypatch.setattr("builtins.input", lambda prompt="": next(answers))
    console.db_catalog()
    out = capsys.readouterr().out
    assert "Hemoglobina, subunitat beta" in out and "Sabies que" in out
    assert "La variant de l'anèmia falciforme" in out
    assert "Anèmia falciforme (drepanocitosi)" in out


def test_texts_are_packaged():
    """La web instala el núcleo desde GitHub: los JSON de textos deben ir en el paquete."""
    import tomllib
    from pathlib import Path
    config = tomllib.loads((Path(__file__).parents[1] / "pyproject.toml").read_text(encoding="utf-8"))
    assert "texts/*.json" in config["tool"]["setuptools"]["package-data"]["dnainvent.database"]
