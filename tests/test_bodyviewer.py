"""El visor 3D de «Cuerpo humano» (bodyviewer) con el cuerpo de juguete de test_anatomy.py.
Lo que dibuja solo se prueba si hay OpenGL sin pantalla (EGL); si no, se salta."""

import numpy as np
import pytest

from dnainvent import anatomy
from dnainvent.anatomy import bodyparts3d, cards
from dnainvent.anatomy.build import build_anatomy
from test_anatomy import asctb, bp_dir, uberon  # noqa: F401  (fixtures)

pytest.importorskip("PIL")
from bodyviewer.camera import OrbitCamera  # noqa: E402
from bodyviewer.info import selection_info  # noqa: E402
from bodyviewer.scene import load_scene, to_gl  # noqa: E402


@pytest.fixture
def built(bp_dir, asctb, uberon, tmp_path, monkeypatch):  # noqa: F811
    out = tmp_path / "anatomy"
    build_anatomy(bodyparts3d.BodyParts3D(bp_dir), asctb, uberon, out, ratio=1.0)
    monkeypatch.setenv("DNAINVENT_ANATOMY", str(out))
    conn = anatomy.open_anatomy()
    yield conn, load_scene(anatomy.meshes_dir())
    conn.close()


def test_coordinates():
    # BodyParts3D: z arriba, y atrás -> OpenGL: y arriba, z hacia quien mira
    assert to_gl([[1, 2, 3]]).tolist() == [[1, 3, -2]]


def test_scene(built):
    _, scene = built
    assert scene.systems == ["skeletal", "heart", "arteries"] and len(scene.pieces) == 4
    assert {p["id"] for p in scene.pieces} == {"FJ1", "FJ2", "FJ3", "FJ4"}
    heart = scene.meshes[1]
    assert set(heart.piece.astype(int)) == {scene.by_id["FJ1"], scene.by_id["FJ4"]}
    assert heart.indices.max() < len(heart.positions)
    assert np.allclose(np.linalg.norm(heart.normals, axis=1), 1, atol=1e-3)
    lo, hi = scene.bounds([scene.by_id["FJ2"]])
    assert np.allclose(hi - lo, 10, atol=0.05)                 # el cubo de 10 mm


def test_camera_frames_a_box():
    cam = OrbitCamera(900)
    cam.frame((-10, -10, -10), (10, 10, 10))
    assert np.allclose(cam.target, 0) and cam.distance < 100
    view = cam.view()
    assert np.allclose(view @ np.array([*cam.eye, 1]), [0, 0, 0, 1], atol=1e-3)
    cam.orbit(10000, 10000)
    assert cam.pitch <= OrbitCamera.MAX_PITCH
    cam.reset()
    assert cam.distance > 900


def test_selection_info(built):
    from i18n import Translator
    conn, _ = built
    card = cards.concept_card(conn, "FMA7088", "es", series={"MYH6"})
    info = selection_info(card, Translator("es"))
    assert info["title"] == "Corazón" and info["subtitle"] == "Heart"
    text = " ".join(t for t, _ in info["lines"])
    assert "★MYH6" in text and "TNNT2" in text and "Corazón" in text
    assert ("★ Genes de serie: MYH6", "gold") in info["lines"]
    assert selection_info(None, Translator("es")) is None


def gl_context():
    moderngl = pytest.importorskip("moderngl")
    try:
        return moderngl.create_standalone_context(backend="egl")
    except Exception:
        pytest.skip("sin OpenGL sin pantalla (EGL)")


def test_render_and_pick(built):
    from bodyviewer.renderer import HIDDEN, SELECTED, Renderer
    _, scene = built
    ctx = gl_context()
    r = Renderer(ctx, scene)
    cam = OrbitCamera(scene.radius)
    femur = scene.by_id["FJ2"]
    cam.frame(*scene.bounds([femur]))
    w, h = 64, 64
    fbo = ctx.simple_framebuffer((w, h), components=4)
    r.render(cam, w, h, fbo)
    assert len(set(fbo.read(components=3))) > 1                      # algo se ha dibujado
    assert r.pick(cam, w, h, w // 2, h // 2) == femur
    r.toggle_system("skeletal")                                        # apagar la capa del fémur
    assert r.state[femur, 3] == HIDDEN and r.pick(cam, w, h, w // 2, h // 2) != femur
    r.select([femur])
    r.isolate = True
    r.update_state()
    assert r.state[femur, 3] == SELECTED and (r.state[:, 3] == HIDDEN).sum() == 3
    assert r.pick(cam, w, h, w // 2, h // 2) == femur
