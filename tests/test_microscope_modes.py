"""Tests de los modos del microscopio: didáctico, realista, atómico (AtomicWorld) y plástico."""

import numpy as np
import pytest

pytest.importorskip("atomicworld")

from atomicworld import atoms
from atomicworld.scene import orbitals
from dnainvent.scene.elements import atom_type
from dnainvent.scene.molecule import build_molecule
from i18n import Translator, available_languages

LANGUAGES = sorted(available_languages())
MODES = ("didactic", "realistic", "atomic", "plastic")


@pytest.fixture(scope="module")
def mol():
    return build_molecule("ATGCGTACGTTAGC")


# ---------- Textos ----------

@pytest.mark.parametrize("language", LANGUAGES)
def test_every_mode_has_its_texts(language):
    tr = Translator(language)
    keys = [f"microscope.hud.mode_{m}" for m in MODES] + [f"microscope.hud.now_{m}" for m in MODES]
    keys += ["microscope.hud.help_atomic", "microscope.hud.only_atomic", "microscope.hud.only_p",
             "microscope.atom.atomic.subshell_title", "microscope.atom.atomic.orbital_title",
             "microscope.atom.atomic.shape.s", "microscope.atom.atomic.shape.p"]
    for key in keys:
        assert tr(key) != key, key
    for key in ("microscope.helix.explanation_atomic", "microscope.helix.explanation_plastic",
                "microscope.atom.atomic.cloud_explanation", "microscope.atom.plastic.explanation",
                "microscope.nucleus.plastic", "microscope.nucleon.plastic"):
        assert tr.lines(key) != [key], key


def test_command_line_mode():
    from microscope.__main__ import parse_args, start_mode
    assert start_mode(parse_args(["ACGT"])) == "didactic"
    assert start_mode(parse_args(["ACGT", "--realistic"])) == "realistic"
    assert start_mode(parse_args(["ACGT", "--mode", "plastic"])) == "plastic"
    with pytest.raises(SystemExit):
        parse_args(["ACGT", "--mode", "cartoon"])


# ---------- Modo plástico: física de la ilustración (sin GPU) ----------

@pytest.mark.parametrize("symbol", ["H", "C", "N", "O", "P"])
def test_one_orbit_per_electron_at_the_subshell_radius(symbol):
    from microscope.levels import plastic
    t = atom_type(symbol)
    o = plastic.electron_orbits(t.element.protons, t.shells, np.random.default_rng(0))
    assert o.count == t.electrons
    expected = np.concatenate([[orbitals.mean_radius_pm(s.n, s.l, atoms.slater_zeff(
        t.element.protons, t.shells, s.n, s.l))] * s.electrons for s in t.shells])
    assert np.allclose(o.radius, expected, rtol=0.09)
    for t_ in (0.0, 1.3):                       # en su círculo todo el tiempo
        r = np.linalg.norm(o.positions(t_), axis=1)
        assert np.allclose(r, o.radius)
    # los internos giran más deprisa (3.ª ley de Kepler), sin parpadear
    w = np.abs(o.omega)
    assert w[0] >= w[-1] and w.max() <= 2 * np.pi / plastic.INNER_PERIOD_S + 1e-9


def test_particle_images_are_there():
    from PIL import Image
    from microscope.render import sprites
    for path in sprites.texture_paths():
        im = Image.open(path)
        assert im.mode == "RGBA" and im.size == (256, 256)
        assert im.getpixel((0, 0))[3] == 0 and im.getpixel((128, 128))[3] == 255   # bola, fondo transparente


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


def _fbo(gl, w=320, h=180):
    fbo = gl.framebuffer(color_attachments=[gl.texture((w, h), 4)],
                         depth_attachment=gl.depth_renderbuffer((w, h)))
    fbo.use()
    return fbo


def _enter_atom(eng, mol, name="N7"):
    i = next(i for i in range(mol.n_atoms) if mol.names[i] == name and mol.pair_index[i] == 2)
    eng.ship.pos = mol.positions[i] + np.array([0.0, 0.0, 0.02])
    eng.ship.vel[:] = 0
    eng.update(1 / 60, np.zeros(3), False)
    assert eng.level.depth == 1
    return i


def test_tab_cycles_through_the_four_modes(gl, mol):
    from microscope.engine import Engine
    eng = Engine(gl, mol, Translator("es"))
    _fbo(gl)
    seen = [eng.mode]
    for _ in range(4):
        seen.append(eng.toggle_mode())
        eng.render(320, 180)
        assert eng.hud_state(320, 180)["mode"] == eng.mode
    assert seen == ["didactic", "realistic", "atomic", "plastic", "didactic"]
    assert eng.toggle_mode(-1) == "plastic" and eng.didactic
    assert eng.set_mode("realistic") == "realistic" and not eng.didactic
    with pytest.raises(ValueError):
        eng.set_mode("cartoon")


@pytest.mark.parametrize("mode", MODES)
def test_full_journey_in_every_mode(gl, mol, mode):
    """Mayús+Enter desde la hélice hasta un quark, dibujando por el camino."""
    from microscope.engine import Engine
    eng = Engine(gl, mol, Translator("ca"), mode=mode)
    _fbo(gl)
    eng.autopilot.dive(eng, full=True)
    depths = set()
    for frame in range(12000):
        eng.update(1 / 60, np.zeros(3), False)
        depths.add(eng.level.depth)
        if frame % 400 == 0:
            eng.render(320, 180)
            eng.hud_state(320, 180)
        if not eng.autopilot.active:
            break
    assert depths == {0, 1, 2, 3}
    eng.render(320, 180)
    assert eng.hud_state(320, 180)["focus"] is not None


def test_plastic_atom_has_a_ball_per_electron_and_nucleon(gl, mol):
    from microscope.engine import Engine
    eng = Engine(gl, mol, Translator("es"), mode="plastic")
    _fbo(gl)
    _enter_atom(eng, mol)
    lvl = eng.level
    atoms_shown = [lvl.i] + lvl.neighbors
    electrons = sum(atom_type(mol.elements[j]).electrons for j in atoms_shown)
    nucleons = sum(atom_type(mol.elements[j]).isotope.a for j in atoms_shown)
    assert lvl.balls.count == electrons + nucleons
    kinds = lvl.balls.data["k"]
    assert (kinds == 2).sum() == electrons                        # capa 2 = electrón
    before = lvl.balls.positions[:electrons].copy()
    eng.update(0.5, np.zeros(3), False)
    assert not np.allclose(before, lvl.balls.positions[:electrons])   # los electrones se mueven
    focus = lvl.focus(np.array([60.0, 0, 0]), np.array([0.0, 1, 0]), True)
    assert any("Bohr" in line for line in focus.explanation)


def test_plastic_balls_are_drawn(gl):
    """Las bolas salen en pantalla, con su color (el protón, rojo) y tapando a las de detrás."""
    from microscope.levels.nucleus import NucleusLevel
    from microscope.mathutil import perspective, view_matrix
    from microscope.render.points import Camera, Renderer
    from microscope.ship import Ship
    fbo = _fbo(gl, 200, 200)
    gl.viewport = (0, 0, 200, 200)
    r = Renderer(gl)
    lvl = NucleusLevel(atom_type("O"), Translator("es"))
    lvl.mode = "plastic"
    lvl.build(r, True)
    ship = Ship(np.array([0.0, 0.0, 9.0]), look_at=np.zeros(3))
    r.begin()
    lvl.render(r, Camera(view_matrix(ship.pos, ship.orient), perspective(70, 1, 0.1, 1e4), 200, 0.0),
               True, 1.0)
    img = np.frombuffer(fbo.read(components=3), dtype=np.uint8).reshape(200, 200, 3).astype(int)
    lit = img.sum(axis=2) > 60
    assert 0.05 < lit.mean() < 0.9
    red = (img[..., 0] > 180) & (img[..., 1] < 120)
    assert red.sum() > 200                                        # algún protón a la vista
    lvl.release()


def test_atomic_mode_keys_o_and_i(gl, mol):
    from microscope.engine import Engine
    eng = Engine(gl, mol, Translator("es"))
    _fbo(gl)
    _enter_atom(eng, mol)                                          # nitrógeno: 1s² 2s² 2p³
    assert eng.cycle_subshell() is None                            # solo en modo atómico
    eng.set_mode("atomic")
    assert eng.cycle_orbital() is None                             # antes hay que elegir un p
    assert eng.cycle_subshell() == "Subnivel 2p³ (O: siguiente)"
    lvl = eng.level
    labels = [eng.cycle_orbital() for _ in range(4)]
    assert labels == ["Orbital 2px (I: siguiente)", "Orbital 2py (I: siguiente)",
                      "Orbital 2pz (I: siguiente)", "Subnivel 2p³ (O: siguiente)"]
    eng.cycle_orbital()
    focus = lvl.focus(np.array([150.0, 0, 0]), np.array([-1.0, 0, 0]), True)
    assert focus.title == "Orbital 2px · 1 electrón" and "[px ↑]" in focus.detail
    eng.render(320, 180)
    assert [name for _, name in lvl.axis_labels()] == ["x", "y", "z"]
    assert [eng.cycle_subshell() for _ in range(3)] == [
        "Subnivel 2s² (O: siguiente)", "Subnivel 1s² (O: siguiente)", "Todos los electrones"]


def test_atomic_mode_lets_you_back_out_to_see_the_whole_atom(gl, mol):
    from microscope.levels.atom import AtomLevel
    from microscope.render.points import Renderer
    i = next(i for i in range(mol.n_atoms) if mol.names[i] == "N7")
    lvl = AtomLevel(mol, i, Translator("es"))
    assert lvl.exit_radius() == pytest.approx(0.8 * lvl.vdw_pm)  # didáctico: desde dentro
    lvl.mode = "atomic"
    assert lvl.exit_radius() > 2 * lvl.vdw_pm
    lvl.build(Renderer(gl), True)
    # todas las capas visibles: la externa (2s, 2p) con más puntos por electrón que la 1s
    assert lvl.clouds["main"].count > 150_000
    lvl.release()


def test_atomic_helix_uses_shell_colors(gl, mol):
    from microscope.levels.atomic_style import SHELL_COLORS, molecule_cloud
    pos, col, size, seed = molecule_cloud(mol, atom_type, 60_000, np.random.default_rng(1))
    assert len(pos) == len(col) == len(size) == len(seed) >= 60 * mol.n_atoms
    rgb = {tuple(np.round(c, 2)) for c in col[:, :3]}
    assert rgb <= {tuple(np.round(SHELL_COLORS[n], 2)) for n in (1, 2, 3)}
