"""Abre el microscopio 3D en un proceso aparte: la consola sigue funcionando mientras tanto."""
import importlib.util
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
MICROSCOPE_PACKAGES = ("numpy", "moderngl", "pyglet", "PIL", "atomicworld")


def missing_dependencies() -> str:
    """Paquetes del microscopio que no están instalados ('' si no falta ninguno)."""
    return ", ".join(p for p in MICROSCOPE_PACKAGES if importlib.util.find_spec(p) is None)


MAX_ARGUMENT = 30_000      # Windows no admite líneas de comando de más de ~32.000 caracteres


def microscope_command(dna: str, start: int, count: int, language: str,
                       realistic: bool = False) -> list[str]:
    """start desde 0, como en el núcleo; la línea de comandos cuenta desde 1, como en biología.

    Si la secuencia es enorme solo se pasa el tramo (y entonces se numera desde 1)."""
    if len(dna) > MAX_ARGUMENT:
        dna, start = dna[start:start + count], 0
    command = [sys.executable, "-m", "microscope", dna, "--start", str(start + 1),
               "--count", str(count), "--lang", language]
    return command + (["--realistic"] if realistic else [])


def open_microscope(dna: str, start: int, count: int, language: str) -> subprocess.Popen:
    return subprocess.Popen(microscope_command(dna, start, count, language), cwd=PROJECT_ROOT)


def protein_microscope_command(session: str, language: str, realistic: bool = False) -> list[str]:
    command = [sys.executable, "-m", "microscope", "--session", str(session), "--lang", language]
    return command + (["--realistic"] if realistic else [])


def open_protein_microscope(session: str, language: str) -> subprocess.Popen:
    return subprocess.Popen(protein_microscope_command(session, language), cwd=PROJECT_ROOT)
