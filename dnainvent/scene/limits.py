"""Cuántos pares de bases enseña el microscopio de una vez (sin NumPy: lo usa la consola).

Una secuencia más larga no se rechaza: se elige un tramo (desde la base N, tantos pares).
Lo mismo con las proteínas: un tramo de aminoácidos (dnainvent/scene/protein.py).
"""
from ..errors import InvalidParameterError

RECOMMENDED_PAIRS = 60     # escritorio: con más, cada átomo lleva menos puntos y va lento
MAX_WEB_PAIRS = 100        # web: máximo por petición a /api/structure (~75 kB comprimido)
# Proteínas: ~15 átomos por aminoácido con los hidrógenos (el ADN, ~63 por par). 150
# aminoácidos son ~2.300 átomos, y al comparar se dibujan dos proteínas a la vez.
RECOMMENDED_RESIDUES = 150


def window(length: int, start: int = 0, count: int | None = None) -> tuple[int, int]:
    """Tramo [start, end) de una secuencia de 'length' elementos (bases o aminoácidos),
    comprobando los límites."""
    if not 0 <= start < length:
        raise InvalidParameterError("start_out_of_range", start=start, length=length)
    if count is not None and count < 1:
        raise InvalidParameterError("count_too_small", minimum=1)
    end = length if count is None else min(length, start + count)
    return start, end
