"""
Prueba sin ventana del microscopio: hace el viaje completo con el piloto automático
(Mayús+Enter) y guarda capturas PNG de cada nivel.

    python tools/microscope_render_test.py [carpeta_salida] [--lang es|en|ca] [--realistic]
    python tools/microscope_render_test.py test_frames --mode plastic     (o atomic)
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import moderngl
import numpy as np
from PIL import Image

from dnainvent.scene.molecule import build_molecule
from i18n import Translator
from microscope.engine import Engine
from microscope.hud import HUD

W, H = 1280, 720


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("out", nargs="?", default="test_frames")
    parser.add_argument("--lang", default="es")
    parser.add_argument("--realistic", action="store_true")
    parser.add_argument("--mode", choices=("didactic", "realistic", "atomic", "plastic"))
    parser.add_argument("--sequence", default="ATGCGTACGTTAGC")
    args = parser.parse_args(argv)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    try:
        ctx = moderngl.create_standalone_context()
    except Exception:
        ctx = moderngl.create_standalone_context(backend="egl")
    fbo = ctx.framebuffer(color_attachments=[ctx.texture((W, H), 4)],
                          depth_attachment=ctx.depth_renderbuffer((W, H)))   # modo plástico
    tr = Translator(args.lang)
    mol = build_molecule(args.sequence)
    eng = Engine(ctx, mol, tr, mode=args.mode or ("realistic" if args.realistic else "didactic"))
    hud = HUD(ctx, tr)

    def snap(name):
        fbo.use()
        eng.render(W, H)
        hud.upload(hud.compose(eng.hud_state(W, H), W, H, eng.time))
        hud.draw()
        img = Image.frombytes("RGBA", (W, H), fbo.read(components=4)).transpose(
            Image.Transpose.FLIP_TOP_BOTTOM)
        img.convert("RGB").save(out / f"{name}.png")
        st = eng.hud_state(W, H)
        print(f"{name}: {eng.level.name!r} · {st['depth_m']:.2e} m · "
              f"{st['focus'].title if st['focus'] else '-'} · {st['autopilot']}")

    snap("0_helix")
    # apuntar a un nitrógeno de la base del par 3 y lanzar el viaje completo
    target = next(i for i in range(mol.n_atoms)
                  if mol.pair_index[i] == 2 and mol.names[i] == "N7" and mol.strand[i] == 1)
    s = eng.ship
    f = mol.positions[target] - s.pos
    f /= np.linalg.norm(f)
    right = np.cross(f, s.up)
    right /= np.linalg.norm(right)
    s.orient = np.column_stack([right, np.cross(right, f), -f])
    eng.autopilot.dive(eng, full=True)
    print("autopiloto:", eng.autopilot.status())

    dt = 1 / 60
    last_depth, shots = 0, 0
    for frame in range(1, 20000):
        eng.update(dt, np.zeros(3), boost=False)
        if eng.level.depth != last_depth:
            last_depth = eng.level.depth
            eng_frames = frame
            shots += 1
        if shots and frame == eng_frames + 90:
            snap(f"{last_depth}_{type(eng.level).__name__}")
        if not eng.autopilot.active:
            snap(f"{last_depth}_end")
            print("fin:", eng.autopilot.message or "-", "en el frame", frame)
            break
    print("Capturas en", out.resolve())


if __name__ == "__main__":
    main()
