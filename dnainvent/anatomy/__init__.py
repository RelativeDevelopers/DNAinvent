"""«Cuerpo humano»: el cuerpo de referencia en 3D, pieza a pieza, con las células de cada parte y
los genes que las marcan.

    from dnainvent import anatomy
    conn = anatomy.open_anatomy()            # ~/.dnainvent/anatomy/anatomy.db, solo lectura
    anatomy.meshes_dir()                     # ~/.dnainvent/anatomy/meshes/<sistema>.dnab

Piezas:   bodyparts3d.py (BodyParts3D 4.0, CC BY 4.0)
Sistemas: systems.py (capas del visor)
Células:  cells.py (tablas ASCT+B del Human Reference Atlas, CC BY 4.0) y tissues.py (tejido
          general escrito a mano)
Mallas:   meshes.py (formato DNAB, el mismo para la consola y la web)
Se construye con tools/build_anatomy.py (build.py). Es divulgación sobre un cuerpo de
referencia (hombre adulto): no sirve para diagnosticar.
"""

import os
import sqlite3
from pathlib import Path

FORMAT_VERSION = 2          # 2: concept_pieces (19-11-2026)


def anatomy_dir() -> Path:
    """DNAINVENT_ANATOMY o ~/.dnainvent/anatomy."""
    return Path(os.environ.get("DNAINVENT_ANATOMY") or Path.home() / ".dnainvent" / "anatomy")


def meshes_dir() -> Path:
    return anatomy_dir() / "meshes"


def database_file() -> Path:
    return anatomy_dir() / "anatomy.db"


def selection_file() -> Path:
    """Lo último que se ha elegido en el visor 3D de escritorio (la consola lo lee)."""
    return anatomy_dir() / "selection.json"


class AnatomyMissing(Exception):
    """Todavía no se ha construido (tools/build_anatomy)."""


def open_anatomy(folder=None) -> sqlite3.Connection:
    path = Path(folder or anatomy_dir()) / "anatomy.db"
    if not path.is_file():
        raise AnatomyMissing(str(path))
    conn = sqlite3.connect(f"file:{path.as_posix()}?mode=ro", uri=True, check_same_thread=False)
    conn.row_factory = sqlite3.Row
    version = conn.execute("SELECT value FROM meta WHERE key = 'format'").fetchone()
    if not version or int(version[0]) != FORMAT_VERSION:
        conn.close()
        raise AnatomyMissing(str(path))
    return conn
