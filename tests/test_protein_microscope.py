"""Tests de la fase 5a: la proteína en el microscopio (nivel, comparación, sesión y lanzador).

Los textos y la navegación no necesitan tarjeta gráfica. El motor (pantalla partida y
superpuestas) sí: esos tests se saltan si no se puede crear un contexto OpenGL.
"""

import random
import re
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("atomicworld")
pytest.importorskip("moderngl")          # los niveles importan el renderizador

from dnainvent.protein import pipeline
from dnainvent.protein.cache import Cache
from dnainvent.protein.clients.blast import SWISSPROT
from dnainvent.protein.compare import compare, moved
from dnainvent.protein.pdb import Structure, read_pdb
from dnainvent.protein.session import Panel, panels_for, read_session, write_session
from dnainvent.scene.elements import atom_type
from dnainvent.scene.limits import RECOMMENDED_RESIDUES
from dnainvent.scene.protein import build_protein
from i18n import Translator, available_languages
from test_protein_compare import SICKLE, FakePredictor, as_prediction
from test_protein_search import HBB, FakeAlphaFold, FakeBlast, blast_hit, gene
from ui import launcher

DATA = Path(__file__).parent / "data"
LANGUAGES = sorted(available_languages())
PLACEHOLDER = re.compile(r"\{(\w+)\}")


def clean(text: str) -> bool:
    return not PLACEHOLDER.search(text) and not re.search(r"\b(microscope|protein|units)\.[a-z_]+", text)


def texts_of(focus) -> list[str]:
    return [focus.title, focus.detail, *focus.explanation] if focus else []


@pytest.fixture(scope="module")
def natural():
    return read_pdb(DATA / "AF-P68871-F1-model_v6.pdb")


@pytest.fixture(scope="module")
def pair(natural):
    """La natural y "la tuya" (Glu7Val, girada y con el extremo C desplazado), ya
    superpuestas y con el mismo origen, como las prepara la sesión."""
    pred = as_prediction(natural, SICKLE, random.Random(5))
    pred = Structure([replace(a, x=a.x + (3.0 if a.resi > 135 else 0.0)) for a in pred.atoms],
                     pred.title)
    c = compare(SICKLE, HBB, natural_structure=natural, predicted=pred)
    origin = tuple(np.mean([a.xyz for a in natural.atoms], axis=0))
    return (build_protein(natural, origin=origin),
            build_protein(moved(pred, c.superposition), origin=origin), c)


# ---------- Nivel de proteína (sin GPU) ----------

@pytest.mark.parametrize("language", LANGUAGES)
def test_protein_level_texts_in_every_language(language, pair):
    from microscope.levels.atom import AtomLevel
    from microscope.levels.nucleus import NucleusLevel
    from microscope.levels.protein import OverlayLevel, ProteinLevel

    tr = Translator(language)
    mol, other, _ = pair
    level = ProteinLevel(mol, tr, "Hemoglobin subunit beta · P68871", {6: "Glu7Val"})
    far = np.array([40.0, 0, 0])
    texts = [level.name, level.time_info()]
    texts += texts_of(level.focus(far, -far / 40, True))
    texts += texts_of(level.focus(far, -far / 40, False))
    for name in ("CA", "CB", "SD", "HA"):
        i = next(k for k, n in enumerate(mol.names) if n == name)
        texts += texts_of(level.atom_focus(i, 1.0, True))
        texts += texts_of(level.atom_focus(i, 1.0, False))
        texts.append(level.dive_target(mol.positions[i] + [0, 0, 1], np.array([0, 0, -1.0]),
                                       True).label())
    mutated = level.ca[6]
    assert any("Glu7Val" in t for t in texts_of(level.atom_focus(mutated, 1.0, True)))

    sulfur = mol.elements.index("S")
    atom = AtomLevel(mol, sulfur, tr)
    texts += [atom.name, *texts_of(atom.focus(np.array([50.0, 0, 0]), np.array([0, 1.0, 0]), True))]
    nucleus = NucleusLevel(atom_type("S"), tr, context=lambda: "x")
    texts += [nucleus.name, *texts_of(nucleus.focus(np.array([20.0, 0, 0]), np.array([0, 1.0, 0]), True))]

    overlay = OverlayLevel([ProteinLevel(mol, tr, "natural", {6: "Glu7Val"}),
                            ProteinLevel(other, tr, "tuya", {6: "Glu7Val"})], tr)
    texts += [overlay.name, *texts_of(overlay.focus(far, -far / 40, True))]
    near = mol.positions[level.ca[20]] + [0.4, 0, 0]
    texts += texts_of(overlay.focus(near, np.array([-1.0, 0, 0]), True))

    bad = [t for t in texts if not clean(t)]
    assert not bad, bad


def test_atom_context_of_a_protein(pair):
    from microscope.levels.helix import atom_context
    mol = pair[0]
    tr = Translator("es")
    i = next(k for k, n in enumerate(mol.names) if n == "CA" and mol.pair_index[k] == 6)
    assert atom_context(tr, mol, i) == "esqueleto · Glutamato (Glu) 7"
    j = next(k for k, n in enumerate(mol.names) if n == "SD")
    assert atom_context(tr, mol, j).startswith("cadena lateral · Metionina (Met)")


def test_dna_atom_context_is_unchanged():
    from dnainvent.scene.molecule import build_molecule
    from microscope.levels.helix import atom_context
    mol = build_molecule("ATGC")
    assert "par 1" in atom_context(Translator("es"), mol, 0)


def test_overview_from_afar_and_atoms_up_close(pair):
    from microscope.levels.protein import ProteinLevel
    mol = pair[0]
    level = ProteinLevel(mol, Translator("en"), "HBB")
    assert level.focus(np.array([30.0, 0, 0]), np.array([-1.0, 0, 0]), True).title == "Protein · HBB"
    i = level.ca[30]
    focus = level.focus(mol.positions[i] + [0.5, 0, 0], np.array([-1.0, 0, 0]), True)
    assert focus.pair_index is not None


def test_entering_an_atom(pair):
    from microscope.levels.protein import OverlayLevel, ProteinLevel
    mol, other, _ = pair
    tr = Translator("es")
    level = ProteinLevel(mol, tr)
    i = mol.elements.index("S")
    child = level.child_to_enter(mol.positions[i] + [0.01, 0, 0])
    assert child.key == ("atom", i) and child.scale == 1000.0
    overlay = OverlayLevel([ProteinLevel(mol, tr), ProteinLevel(other, tr)], tr)
    j = other.elements.index("S")
    child = overlay.child_to_enter(other.positions[j] + [0.001, 0, 0])
    assert child is not None and child.key[-1] in (0, 1)
    target = overlay.dive_target(other.positions[j] + [0.5, 0, 0], np.array([-1.0, 0, 0]), True)
    assert len(target.key) == 3


def test_both_proteins_share_coordinates(pair):
    mol, other, c = pair
    ca_mol = {int(mol.pair_index[i]): mol.positions[i] for i, n in enumerate(mol.names) if n == "CA"}
    ca_other = {int(other.pair_index[i]): other.positions[i]
                for i, n in enumerate(other.names) if n == "CA"}
    d = np.array([np.linalg.norm(ca_mol[k] - ca_other[k]) for k in range(130)]) * 10   # Å
    assert d.max() < 1.0                         # superpuestas: la misma forma, el mismo sitio


def test_build_protein_with_a_given_origin(natural):
    origin = (0.0, 0.0, 0.0)
    mol = build_protein(natural, origin=origin)
    ca = next(a for a in natural.residues()[0].atoms.values() if a.name == "N")
    assert np.allclose(mol.positions[0], np.array(ca.xyz) / 10, atol=1e-4)


# ---------- Sesión ----------

def sickle_report(natural, with_structure=True, with_prediction=True):
    hit = replace(blast_hit(identity=99.3), query_aligned=SICKLE, hit_aligned=HBB)
    afdb = FakeAlphaFold() if with_structure else FakeAlphaFold(
        __import__("dnainvent").ServiceError("protein_not_found", accession="P68871"))
    predictor = FakePredictor(natural) if with_prediction else None
    return pipeline.analyze(gene(SICKLE), FakeBlast({SWISSPROT: [hit]}), afdb, predictor)


LABELS = {"natural": "HBB · P68871", "yours": "Tu proteína (ESMFold)"}


def test_session_for_a_variant(tmp_path, natural):
    report = sickle_report(natural)
    panels = panels_for(report, LABELS)
    assert [p.label for p in panels] == list(LABELS.values())
    assert panels[0].highlight == {6: "Glu7Val"} and panels[1].highlight == {6: "Glu7Val"}
    path = write_session(panels, Cache(tmp_path))
    data = read_session(path)
    assert len(data["panels"]) == 2
    assert data["panels"][0]["highlight"] == {6: "Glu7Val"}
    assert data["panels"][1]["structure"].sequence() == SICKLE
    assert len(data["origin"]) == 3


def test_session_puts_yours_on_top_of_the_natural(tmp_path, natural):
    report = sickle_report(natural)
    data = read_session(write_session(panels_for(report, LABELS), Cache(tmp_path)))
    first, second = (p["structure"].residues() for p in data["panels"])
    gap = max(np.linalg.norm(np.subtract(a.atoms["CA"].xyz, b.atoms["CA"].xyz))
              for a, b in zip(first, second))
    assert gap < 0.01                           # la predicción falsa es la natural girada


def test_session_known_protein_has_one_panel(natural):
    report = pipeline.analyze(gene(HBB), FakeBlast({SWISSPROT: [blast_hit()]}), FakeAlphaFold(),
                              FakePredictor(natural))
    panels = panels_for(report, LABELS)
    assert [p.label for p in panels] == ["HBB · P68871"] and panels[0].highlight == {}


def test_session_invented_protein_shows_yours(natural):
    report = pipeline.analyze(gene(HBB), FakeBlast({}), None, FakePredictor(natural))
    assert [p.label for p in panels_for(report, LABELS)] == ["Tu proteína (ESMFold)"]


def test_session_without_natural_structure(natural):
    report = sickle_report(natural, with_structure=False)
    panels = panels_for(report, LABELS)
    assert [p.label for p in panels] == ["Tu proteína (ESMFold)"]
    assert panels[0].highlight == {6: "Glu7Val"}


def test_session_window_maps_to_the_other_protein(natural):
    report = sickle_report(natural)
    panels = panels_for(report, LABELS, start=100, count=30)
    assert (panels[0].start, panels[0].count) == (100, 30)
    assert (panels[1].start, panels[1].count) == (100, 30)


def test_microscope_loads_a_session(tmp_path, natural):
    from microscope.__main__ import load_protein, parse_args
    path = write_session(panels_for(sickle_report(natural), LABELS), Cache(tmp_path))
    mol, comparison, label, highlight = load_protein(parse_args(["--session", str(path)]),
                                                     Translator("es"))
    assert mol.kind == "protein" and comparison.second.sequence == SICKLE
    assert comparison.labels == tuple(LABELS.values()) and highlight == {6: "Glu7Val"}


def test_microscope_loads_a_pdb():
    from microscope.__main__ import load_protein, parse_args
    args = parse_args(["--pdb", str(DATA / "AF-P68871-F1-model_v6.pdb"), "--start", "11",
                       "--count", "20"])
    mol, comparison, _, _ = load_protein(args, Translator("es"))
    assert comparison is None and (mol.start, mol.end) == (10, 30)


def test_protein_launcher_command():
    cmd = launcher.protein_microscope_command("/tmp/s.json", "ca", realistic=True)
    assert cmd[:3] == [sys.executable, "-m", "microscope"]
    assert cmd[cmd.index("--session") + 1] == "/tmp/s.json"
    assert cmd[cmd.index("--lang") + 1] == "ca" and "--realistic" in cmd


def test_recommended_residues_is_reasonable():
    assert 50 <= RECOMMENDED_RESIDUES <= 400


# ---------- Motor (necesita OpenGL) ----------

@pytest.fixture(scope="module")
def gl():
    import moderngl
    for kwargs in ({}, {"backend": "egl"}):
        try:
            return moderngl.create_standalone_context(**kwargs)
        except Exception:
            continue
    pytest.skip("sin contexto OpenGL")


def test_engine_split_and_overlay(gl, pair):
    from microscope.engine import Comparison, Engine
    mol, other, _ = pair
    tr = Translator("es")
    eng = Engine(gl, mol, tr, True, Comparison(other, ("natural", "tuya"),
                                               ({6: "Glu7Val"}, {6: "Glu7Val"})))
    fbo = gl.framebuffer(color_attachments=[gl.texture((320, 180), 4)])
    fbo.use()
    assert eng.layout == "split" and eng.split_visible
    eng.render(320, 180)
    state = eng.hud_state(320, 180)
    assert [p["label"] for p in state["panels"]] == ["natural", "tuya"]
    assert state["comparing"] and state["molecule"] is mol

    assert eng.switch_side() and eng.view.label == "tuya"
    assert eng.hud_state(320, 180)["molecule"] is other

    assert eng.toggle_layout() and eng.layout == "overlay"
    eng.render(320, 180)
    assert eng.hud_state(320, 180)["panels"] == []
    assert eng.toggle_layout() and eng.layout == "split"

    # dentro de un átomo: pantalla completa del lado activo, y no se cambia la disposición
    i = other.elements.index("S")
    eng.ship.pos = other.positions[i] + np.array([0.01, 0, 0])
    eng.update(1 / 60, np.zeros(3), False)
    assert eng.level.depth == 1 and not eng.split_visible
    assert not eng.toggle_layout() and not eng.switch_side()
    eng.render(320, 180)


def test_engine_single_protein(gl, pair):
    from microscope.engine import Engine
    eng = Engine(gl, pair[0], Translator("en"), False, label="HBB")
    assert eng.layout == "single" and not eng.toggle_layout()
    fbo = gl.framebuffer(color_attachments=[gl.texture((160, 90), 4)])
    fbo.use()
    eng.render(160, 90)
    assert eng.hud_state(160, 90)["focus"].title == "Protein · HBB"
