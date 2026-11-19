"""Abre el cuerpo humano en 3D.

    python -m bodyviewer                      el cuerpo entero
    python -m bodyviewer --select FMA7088     con el corazón elegido y centrado
    python -m bodyviewer --lang ca

Necesita numpy, moderngl, pyglet y pillow, y que «Cuerpo humano» esté preparado
(tools/build_anatomy; ver README).
"""
import argparse
import sys

from i18n import DEFAULT_LANGUAGE, Translator, available_languages
from ui import settings


def parse_args(argv=None):
    parser = argparse.ArgumentParser(prog="python -m bodyviewer",
                                     description="DNAinvent — cuerpo humano en 3D / human body in 3D")
    parser.add_argument("--select", help="parte del cuerpo elegida al abrir (ID de la FMA, p. ej. FMA7088)")
    parser.add_argument("--lang", choices=sorted(available_languages()), help="es, en o ca")
    return parser.parse_args(argv)


def series_genes() -> frozenset:
    """Los genes de serie de «Genes reales», para marcarlos con ★."""
    try:
        from dnainvent.database.catalog import CATALOG_GENES
        return frozenset(CATALOG_GENES)
    except ImportError:
        return frozenset()


def main(argv=None):
    args = parse_args(argv)
    saved = settings.language()
    tr = Translator(args.lang or (saved if saved in available_languages() else DEFAULT_LANGUAGE))
    from dnainvent import anatomy
    try:
        conn = anatomy.open_anatomy()
    except anatomy.AnatomyMissing as exc:
        sys.exit(tr("anatomy.missing", path=exc))
    from bodyviewer.scene import load_scene
    scene = load_scene(anatomy.meshes_dir())
    import pyglet
    from bodyviewer.app import App
    App(scene, conn, tr, series_genes(), args.select)
    pyglet.app.run()


if __name__ == "__main__":
    main()
