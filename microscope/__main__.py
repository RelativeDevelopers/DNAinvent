"""Arranca el microscopio 3D: un viaje desde la doble hélice (o una proteína) hasta los quarks.

    python -m microscope                              (pide la secuencia)
    python -m microscope ATGCGTACGTTAGC
    python -m microscope ATGC --lang ca --realistic
    python -m microscope <secuencia larga> --start 120 --count 60   (un tramo)

    python -m microscope --pdb proteina.pdb           una proteína (cualquier archivo PDB)
    python -m microscope --pdb grande.pdb --start 200 --count 150   (un tramo de aminoácidos)
    python -m microscope --session sesion.json        lo que prepara la opción 8 de la consola
                                                      (la natural y la tuya, para compararlas)

Necesita numpy, moderngl, pyglet, pillow y AtomicWorld (ver README).
"""
import argparse
import sys

from dnainvent import DNAinventError
from dnainvent.scene.limits import RECOMMENDED_PAIRS, RECOMMENDED_RESIDUES
from dnainvent.validation import validate_dna
from i18n import DEFAULT_LANGUAGE, Translator, available_languages
from ui import settings

EXAMPLE = "ATGCGTACGTTAGC"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(prog="python -m microscope",
                                     description="DNAinvent — microscopio 3D / 3D microscope")
    parser.add_argument("sequence", nargs="?", help="ADN (A, C, G, T), p. ej. ATGCGTACGTTAGC")
    parser.add_argument("--start", type=int, default=1,
                        help="primera base del tramo, desde 1 (secuencias largas)")
    parser.add_argument("--count", type=int, default=None,
                        help=f"pares del tramo (por defecto, hasta {RECOMMENDED_PAIRS})")
    parser.add_argument("--lang", choices=sorted(available_languages()),
                        help="es, en o ca (por defecto, el idioma guardado)")
    parser.add_argument("--realistic", action="store_true", help="empezar en modo realista")
    parser.add_argument("--pdb", help="archivo PDB de una proteína")
    parser.add_argument("--session", help="sesión preparada por la consola (opción 8)")
    return parser.parse_args(argv)


def load_protein(args, tr):
    """(molécula, Comparison o None, nombre, resaltado) para --pdb o --session."""
    from dnainvent.protein.pdb import read_pdb
    from dnainvent.protein.session import read_session
    from dnainvent.scene.protein import build_protein
    from microscope.engine import Comparison

    if args.pdb:
        structure = read_pdb(args.pdb)
        count = args.count or RECOMMENDED_RESIDUES
        molecule = build_protein(structure, args.start - 1, count)
        return molecule, None, structure.title.title()[:60], {}
    session = read_session(args.session)
    origin = tuple(session["origin"])
    built = [build_protein(p["structure"], p["start"], p["count"], origin=origin)
             for p in session["panels"]]
    panels = session["panels"]
    if len(built) == 1:
        return built[0], None, panels[0]["label"], panels[0]["highlight"]
    comparison = Comparison(built[1], (panels[0]["label"], panels[1]["label"]),
                            (panels[0]["highlight"], panels[1]["highlight"]))
    return built[0], comparison, panels[0]["label"], panels[0]["highlight"]


def ask_sequence(tr) -> str:
    while True:
        text = input(tr("microscope.console.ask", example=EXAMPLE)).strip()
        if not text:
            return EXAMPLE
        try:
            return validate_dna(text)
        except DNAinventError as exc:
            print("  " + tr.error(exc))


def main(argv=None):
    args = parse_args(argv)
    saved = settings.language()
    tr = Translator(args.lang or (saved if saved in available_languages() else DEFAULT_LANGUAGE))
    if args.pdb or args.session:
        try:
            molecule, comparison, label, highlight = load_protein(args, tr)
        except (DNAinventError, OSError, ValueError) as exc:
            sys.exit(tr.error(exc) if isinstance(exc, DNAinventError) else str(exc))
        import pyglet
        from microscope.app import App
        App(molecule, tr, not args.realistic, comparison, label, highlight)
        pyglet.app.run()
        return
    try:
        dna = validate_dna(args.sequence) if args.sequence else ask_sequence(tr)
        start = args.start - 1
        count = args.count or min(RECOMMENDED_PAIRS, len(dna) - max(start, 0))
        from dnainvent.scene.molecule import build_molecule
        molecule = build_molecule(dna, start, count)
    except DNAinventError as exc:
        sys.exit(tr.error(exc))
    if count > RECOMMENDED_PAIRS:
        print(tr("microscope.console.slow", pairs=count, maximum=RECOMMENDED_PAIRS))

    import pyglet                       # solo aquí: los tests no necesitan ventana
    from microscope.app import App
    App(molecule, tr, didactic=not args.realistic)
    pyglet.app.run()


if __name__ == "__main__":
    main()
