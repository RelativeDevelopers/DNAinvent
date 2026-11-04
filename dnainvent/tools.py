"""Herramientas del submenú '*': buscar, similitud y permutaciones."""

import math
from collections import Counter
from collections.abc import Iterator

from .alignment import align
from .errors import InvalidParameterError
from .validation import validate_dna


def find_positions(dna: str, target: str) -> list[int]:
    """Posiciones (desde 0) donde aparece target dentro de dna, sin solapamiento."""
    dna, target = validate_dna(dna), validate_dna(target)
    positions, start = [], dna.find(target)
    while start != -1:
        positions.append(start)
        start = dna.find(target, start + len(target))
    return positions


def similarity(first: str, second: str) -> float:
    """Porcentaje de identidad tras alinear las dos secuencias."""
    return align(first, second).identity


def count_positional_permutations(dna: str, r: int) -> int:
    """Ordenaciones de r posiciones, tratando cada base como distinta aunque se repita.
    Es lo que cuenta itertools.permutations: n! / (n - r)!"""
    _check_r(r)
    return math.perm(len(validate_dna(dna)), r)


def count_permutations(dna: str, r: int) -> int:
    """Cuántas secuencias DISTINTAS de r bases se pueden formar, sin generarlas.

    Se añade una base cada vez. Si ya hay 'ways[j]' secuencias de longitud j,
    añadir k copias de la base nueva da secuencias de longitud j + k, y hay
    C(j + k, k) formas de elegir en qué posiciones van esas k copias.
    Ejemplo: "AAC", r = 3 -> 3 (AAC, ACA, CAA), mientras que 3! = 6.
    """
    _check_r(r)
    counts = Counter(validate_dna(dna))
    ways = [1] + [0] * r  # ways[j]: secuencias distintas de longitud j
    for copies in counts.values():
        new_ways = [0] * (r + 1)
        for j, current in enumerate(ways):
            if current:
                for k in range(min(copies, r - j) + 1):
                    new_ways[j + k] += current * math.comb(j + k, k)
        ways = new_ways
    return ways[r]


def iter_permutations(dna: str, r: int) -> Iterator[str]:
    """Genera las permutaciones distintas de r bases, en orden alfabético.

    En lugar de ordenar posiciones (que repite resultados cuando hay bases
    iguales), en cada paso elige qué base poner entre las que aún quedan.
    Como cada base se prueba una sola vez por paso, no salen duplicados.

    Es un backtracking con una pila propia en lugar de recursión: una función
    que se llama a sí misma falla (RecursionError) hacia los 1.000 niveles,
    y aquí cada base de la permutación sería un nivel.
    stack[d] guarda qué base (índice en 'bases') toca probar a continuación
    en la posición d.
    """
    _check_r(r)
    remaining = Counter(validate_dna(dna))
    bases = sorted(remaining)
    if r > sum(remaining.values()):
        return

    current: list[str] = []
    stack = [0]

    def step_back() -> None:
        """Abandona la posición actual y devuelve al contador la base que llevó hasta ella."""
        stack.pop()
        if current:
            remaining[current.pop()] += 1

    while stack:
        if len(current) == r:                  # permutación completa
            yield "".join(current)
            step_back()
            continue
        i = stack[-1]
        while i < len(bases) and not remaining[bases[i]]:
            i += 1                             # saltar bases agotadas
        if i == len(bases):                    # no quedan opciones en esta posición
            step_back()
            continue
        stack[-1] = i + 1                      # la próxima vez, probar la siguiente base
        remaining[bases[i]] -= 1
        current.append(bases[i])
        stack.append(0)                        # bajar a la siguiente posición


def _check_r(r: int) -> None:
    if r < 0:
        raise InvalidParameterError("negative_permutation_size")
