"""Tests del microscopio 3D sin ventana ni tarjeta gráfica: textos, navegación y lanzador."""

import re
import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("atomicworld")

from dnainvent.scene.molecule import build_molecule
from i18n import Translator, available_languages
from ui import launcher, settings

ROOT = Path(__file__).parents[1]
LANGUAGES = sorted(available_languages())
PLACEHOLDER = re.compile(r"\{(\w+)\}")
KEY = re.compile(r'"((?:microscope|units)\.[a-z_.]+[a-z_])"')


@pytest.fixture(scope="module")
def mol():
    return build_molecule("ATGCGTACGTTAGC")


def clean(text: str) -> bool:
    """Un texto bien traducido: ni claves sin traducir ni {huecos} sin rellenar."""
    return not PLACEHOLDER.search(text) and not re.search(r"\b(microscope|units)\.[a-z_]+", text)


def all_texts(focus) -> list[str]:
    return [focus.title, focus.detail, *focus.explanation]


def test_every_microscope_key_used_in_code_exists():
    tr = Translator("es")
    used = set()
    for path in (ROOT / "microscope").rglob("*.py"):
        used |= set(KEY.findall(path.read_text(encoding="utf-8")))
    used |= set(KEY.findall((ROOT / "ui" / "console.py").read_text(encoding="utf-8")))
    missing = {k for k in used if tr(k) == k and tr.lines(k) == [k]}
    # las claves construidas con f-strings (microscope.nucleus.{kind}...) se prueban abajo
    assert not missing, missing


@pytest.mark.parametrize("language", LANGUAGES)
def test_level_texts_in_every_language(language, mol):
    pytest.importorskip("moderngl")    # los niveles importan el renderizador (GPU)
    from microscope.levels.atom import AtomLevel
    from microscope.levels.helix import HelixLevel
    from microscope.levels.nucleon import NucleonLevel
    from microscope.levels.nucleus import NucleusLevel
    from dnainvent.scene.elements import atom_type

    tr = Translator(language)
    helix = HelixLevel(mol, tr)
    far = np.array([30.0, 0.0, 0.0])
    texts = [helix.name, helix.time_info(), *all_texts(helix.focus(far, -far / 30, True))]
    for part in ("base", "sugar", "phosphate"):
        i = mol.part.index(part)
        texts += all_texts(helix.atom_focus(i, 1.0, True))
        texts.append(helix.dive_target(mol.positions[i] + [0, 0, 1], np.array([0, 0, -1.0]), True).label())

    n7 = mol.names.index("N7")
    atom = AtomLevel(mol, n7, tr)
    texts += [atom.name, atom.time_info()]
    texts += all_texts(atom.focus(np.array([50.0, 0, 0]), np.array([-1.0, 0, 0]), True))   # núcleo
    texts += all_texts(atom.focus(np.array([50.0, 0, 0]), np.array([0, 1.0, 0]), True))    # nube
    texts += all_texts(atom.focus(np.array([50.0, 0, 0]), np.array([0, 1.0, 0]), False))
    texts.append(atom.dive_target(np.zeros(3), np.array([1.0, 0, 0]), True).label())

    for symbol in ("H", "C", "P"):
        nucleus = NucleusLevel(atom_type(symbol), tr, context=lambda: "x")
        texts += [nucleus.name, nucleus.time_info()]
        texts += all_texts(nucleus.focus(np.array([20.0, 0, 0]), np.array([0, 1.0, 0]), True))
        texts += all_texts(nucleus.focus(np.array([20.0, 0, 0]), np.array([0, 1.0, 0]), False))
        k = 0
        texts += all_texts(nucleus.focus(nucleus.centers[k] + [0.3, 0, 0], np.array([-1.0, 0, 0]), True))
        texts.append(nucleus.dive_target(np.array([20.0, 0, 0]), np.array([-1.0, 0, 0]), True).label())

    for kind in ("p", "n"):
        nucleon = NucleonLevel(kind, 3, tr)
        texts += [nucleon.name, nucleon.time_info()]
        texts += all_texts(nucleon.focus(np.array([9.0, 0, 0]), np.array([0, 1.0, 0]), True))
        texts += all_texts(nucleon.focus(np.array([9.0, 0, 0]), np.array([0, 1.0, 0]), False))
        q = nucleon.quark_pos[0]
        texts += all_texts(nucleon.focus(q + [0.5, 0, 0], np.array([-1.0, 0, 0]), True))
        texts.append(nucleon.dive_target(np.zeros(3), np.array([1.0, 0, 0]), True).label())

    bad = [t for t in texts if not clean(t)]
    assert not bad, bad


def test_catalan_and_spanish_prepositions():
    assert Translator("es").of("Oxígeno") == "de oxígeno"
    assert Translator("ca").of("Oxigen") == "d'oxigen"
    assert Translator("ca").of("Hidrogen") == "d'hidrogen"
    assert Translator("ca").of("Carboni") == "de carboni"
    assert Translator("en").of("Carbon") == "Carbon"


def test_element_names_come_from_atomicworld():
    assert Translator("es").element(15) == "Fósforo"
    assert Translator("ca").element(15) == "Fòsfor"
    assert Translator("en").element(15) == "Phosphorus"


def test_numbers_and_plurals():
    tr = Translator("ca")
    assert tr.count(1, "proton") == "1 protó"
    assert tr.count(15, "proton") == "15 protons"
    assert tr.sci(2.3e17) == "2,3×10¹⁷"
    assert Translator("en").sci(1.2e-10) == "1.2×10⁻¹⁰"


def test_next_language_cycles_through_all():
    tr = Translator("es")
    seen = {tr.next_language() for _ in LANGUAGES}
    assert seen == set(LANGUAGES)


def test_settings_roundtrip(tmp_path):
    path = tmp_path / "settings.json"
    assert settings.language(path) is None
    settings.save({"language": "ca"}, path)
    assert settings.language(path) == "ca"
    path.write_text("{roto", encoding="utf-8")
    assert settings.language(path) is None          # un archivo roto no impide arrancar


def test_launcher_command_counts_from_one():
    cmd = launcher.microscope_command("ACGT" * 50, start=119, count=60, language="ca")
    assert cmd[:3] == [sys.executable, "-m", "microscope"]
    assert cmd[cmd.index("--start") + 1] == "120"
    assert cmd[cmd.index("--count") + 1] == "60"
    assert cmd[cmd.index("--lang") + 1] == "ca"


def test_launcher_passes_only_the_window_of_huge_sequences():
    dna = "A" * (launcher.MAX_ARGUMENT + 10)
    cmd = launcher.microscope_command(dna, start=500, count=60, language="es")
    assert len(cmd[3]) == 60 and cmd[cmd.index("--start") + 1] == "1"


def test_command_line_arguments():
    from microscope.__main__ import parse_args
    args = parse_args(["ACGT", "--start", "3", "--count", "2", "--lang", "en", "--realistic"])
    assert (args.sequence, args.start, args.count, args.lang, args.realistic) == ("ACGT", 3, 2, "en", True)
