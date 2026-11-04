"""Punto de entrada de DNAinvent.

    python main.py                 -> consola (la primera vez pregunta el idioma y lo guarda)
    python main.py --lang ca       -> consola en catalán (también: es, en)
    python -m microscope ATGC      -> microscopio 3D directamente (ver microscope/__main__.py)

La versión web está en su propio repositorio: ../dnainvent_web
"""

import argparse
import os

from i18n import available_languages


def main() -> None:
    parser = argparse.ArgumentParser(description="DNAinvent")
    parser.add_argument("--lang", choices=sorted(available_languages()),
                        default=os.environ.get("DNAINVENT_LANG"),
                        help="es, en o ca (por defecto, el idioma guardado)")
    args = parser.parse_args()

    from ui.console import run
    run(args.lang)


if __name__ == "__main__":
    main()
