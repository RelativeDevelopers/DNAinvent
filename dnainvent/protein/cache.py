"""Caché en disco: lo descargado o calculado una vez no se vuelve a pedir.

Las búsquedas BLAST tardan, ESMFold limita cuántas peticiones se le hacen y las estructuras
de AlphaFold DB no cambian dentro de una versión. Todo se guarda en ~/.dnainvent/cache
(la misma carpeta de los ajustes), en subcarpetas por servicio:

    ~/.dnainvent/cache/afdb/P68871.json                 ficha de la entrada (caduca)
    ~/.dnainvent/cache/afdb/AF-P68871-F1-v6.pdb          estructura (no caduca)

Para usar otra carpeta: variable de entorno DNAINVENT_CACHE. Borrar la carpeta es seguro.
"""

import hashlib
import os
import re
import time
from pathlib import Path

_SAFE = re.compile(r"[^A-Za-z0-9._-]")


def _safe(name: str) -> str:
    """Nombre de archivo sin barras ni caracteres raros; "." y ".." no se permiten (saldrían
    de la carpeta de la caché)."""
    name = _SAFE.sub("_", name)
    return name if name.strip(".") else "_"


def default_root() -> Path:
    return Path(os.environ.get("DNAINVENT_CACHE") or Path.home() / ".dnainvent" / "cache")


def key_for(text: str) -> str:
    """Clave corta y segura para textos largos (p. ej. una secuencia de 1.000 aminoácidos)."""
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:32]


class Cache:
    def __init__(self, root: str | Path | None = None):
        self.root = Path(root) if root is not None else default_root()

    def path(self, namespace: str, key: str) -> Path:
        return self.root / _safe(namespace) / _safe(key)

    def get(self, namespace: str, key: str, max_age: float | None = None) -> bytes | None:
        """Contenido guardado, o None si no existe o tiene más de max_age segundos."""
        path = self.path(namespace, key)
        try:
            if max_age is not None and time.time() - path.stat().st_mtime > max_age:
                return None
            return path.read_bytes()
        except OSError:
            return None

    def put(self, namespace: str, key: str, data: bytes) -> Path:
        """Guarda de forma atómica: se escribe en un temporal y se renombra, así un corte a
        medias nunca deja un archivo roto en la caché."""
        path = self.path(namespace, key)
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(f"{path.name}.{os.getpid()}.tmp")
        temporary.write_bytes(data)
        os.replace(temporary, path)
        return path
