"""Fichas de «Cuerpo humano» (dnainvent.anatomy.cards y names) y la consola, con el cuerpo de
juguete de test_anatomy.py (un corazón con dos piezas, un fémur y una aorta)."""

import json

import pytest

from dnainvent import anatomy
from dnainvent.anatomy import bodyparts3d, cards, names
from dnainvent.anatomy.build import build_anatomy
from test_anatomy import asctb, bp_dir, uberon  # noqa: F401  (fixtures)


@pytest.fixture
def built(bp_dir, asctb, uberon, tmp_path, monkeypatch):  # noqa: F811
    out = tmp_path / "anatomy"
    build_anatomy(bodyparts3d.BodyParts3D(bp_dir), asctb, uberon, out, ratio=1.0)
    monkeypatch.setenv("DNAINVENT_ANATOMY", str(out))
    conn = anatomy.open_anatomy()
    yield conn
    conn.close()


def test_names():
    assert names.concept("FMA7088", "es") == "Corazón" and names.concept("FMA7088", "ca") == "Cor"
    assert names.concept("FMA7088", "en") is None
    assert names.normalize("  Fémur Cèl·lula ") == "femur cellula"
    assert names.english("left femur") == "Left femur"


def test_search(built):
    hits = cards.search(built, "corazon", "es")              # sin acento también
    assert hits[0]["fma"] == "FMA7088" and hits[0]["label"] == "Corazón" and hits[0]["pieces"] == 2
    assert [h["fma"] for h in cards.search(built, "heart", "en")][0] == "FMA7088"
    assert cards.search(built, "femur", "ca")[0]["name"] == "Left femur"
    assert cards.search(built, "x", "es") == [] and cards.search(built, "nada parecido", "es") == []


def test_heart_card_joins_the_cells_of_its_parts(built):
    card = cards.concept_card(built, "FMA7088", "es", series={"MYH6"})
    assert (card["label"], card["system"], card["pieces"]) == ("Corazón", "heart", 2)
    assert card["cells"]["source"] == "own" and card["cells"]["in_parts"]
    assert {g["symbol"] for g in card["genes"]} == {"TNNT2", "MYH6", "MYH11"}
    assert card["series_genes"] == ["MYH6"]
    assert card["cells"]["cells"][0]["label"] == "célula del músculo cardíaco"       # la que más genes tiene
    assert ("Uberon UBERON:0000948", cards.UBERON_URL.format(id="UBERON:0000948")) in card["links"]
    assert [p["fma"] for p in card["parts"]] == ["FMA7100", "FMA9999"]


def test_part_takes_the_cells_of_its_organ(built):
    card = cards.concept_card(built, "FMA7100", "en")
    assert card["cells"]["source"] == "organ" and card["cells"]["found_at"]["fma"] == "FMA7088"
    assert [p["fma"] for p in card["part_of"]] == ["FMA7088"]


def test_bone_gets_the_general_tissue(built):
    card = cards.concept_card(built, "FMA9612", "ca")
    assert card["cells"]["source"] == "tissue" and card["cells"]["tissue"] == "bone"
    assert "RUNX2" in {g["symbol"] for g in card["genes"]}
    assert ("Cell Ontology", "https://obofoundry.org/ontology/cl.html") in card["links"]
    assert card["size_cm"] == 1.0 and cards.concept_card(built, "FMA0", "es") is None


def test_pieces_and_systems(built):
    assert cards.concept_pieces(built, "FMA7088") == ["FJ1", "FJ4"]
    assert cards.piece_concept(built, "FJ2") == "FMA9612"
    assert [s["key"] for s in cards.systems(built)] == ["skeletal", "heart", "arteries"]
    assert [i["fma"] for i in cards.system_concepts(built, "heart", "es")] == ["FMA7100", "FMA9999"]
    assert cards.system_concepts(built, "heart", "en", "wall")[0]["name"] == "Wall of heart"


def test_size_of_a_pair_is_one_side():
    left = {"lo_x": -50, "hi_x": -10, "lo_y": 0, "hi_y": 10, "lo_z": 0, "hi_z": 400}
    right = {"lo_x": 10, "hi_x": 50, "lo_y": 0, "hi_y": 10, "lo_z": 0, "hi_z": 400}
    assert cards.size_cm([left, right]) == 40.0 and cards.size_cm([]) is None


@pytest.mark.parametrize("language", ["es", "ca", "en"])
def test_console_card(built, language, monkeypatch, capsys):
    from i18n import Translator
    from ui import console
    monkeypatch.setattr(console, "tr", Translator(language))
    monkeypatch.setattr(console, "_anatomy", None)
    answers = iter(["1", "corazon" if language == "es" else "heart", "1", "0", "4", "0", "0"])
    anatomy.selection_file().write_text(json.dumps({"fma": "FMA9612", "piece": "FJ2"}), encoding="utf-8")
    monkeypatch.setattr("builtins.input", lambda prompt="": next(answers))
    console.anatomy_menu()
    out = capsys.readouterr().out
    assert "TNNT2" in out and "RUNX2" in out and "BodyParts3D" in out
    assert "{" not in out.replace("{}", "")                     # ningún texto sin rellenar
    assert {"es": "Corazón", "ca": "Cor", "en": "Heart"}[language] in out


def test_console_without_anatomy(tmp_path, monkeypatch, capsys):
    from i18n import Translator
    from ui import console
    monkeypatch.setenv("DNAINVENT_ANATOMY", str(tmp_path))
    monkeypatch.setattr(console, "tr", Translator("es"))
    monkeypatch.setattr(console, "_anatomy", None)
    console.anatomy_search()
    assert "tools/build_anatomy" in capsys.readouterr().out
