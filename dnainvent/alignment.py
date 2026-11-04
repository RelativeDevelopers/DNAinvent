"""Alineamiento global de dos secuencias (Needleman-Wunsch, 1970; huecos afines de Gotoh, 1982).

Idea: buscar la forma de colocar una secuencia debajo de la otra, insertando
huecos ('-') donde haga falta, para que se parezcan lo máximo posible.

    Sin alinear:   ACGTACGTAC        Alineadas:   -ACGTACGTAC
                   TACGTACGTAC                    TACGTACGTAC

Cada posición del alineamiento suma o resta puntos. Hay dos formas de puntuar:

  ADN (Scoring): coincidencia +1, sustitución -1, hueco -2 por posición.

  Proteínas (MatrixScoring): no todos los cambios valen lo mismo. Cambiar una isoleucina
  por una valina (las dos pequeñas e hidrófobas) es frecuente en la evolución; cambiarla por
  un aspartato (con carga) es raro. La matriz BLOSUM62 da una puntuación a cada pareja.
  Además los huecos son "afines": abrir un hueco cuesta mucho (-12) y alargarlo poco (-1
  por posición), porque en la evolución una inserción de 5 aminoácidos suele ser un único
  suceso, no cinco.

Probar todas las colocaciones posibles sería inabarcable, así que se usa
programación dinámica: se rellena una tabla donde cada casilla guarda la
mejor puntuación posible para alinear los primeros i caracteres de una
secuencia con los primeros j de la otra, reutilizando las casillas ya
calculadas. Coste: proporcional a len(a) * len(b).

Para que los huecos afines funcionen, cada casilla guarda tres puntuaciones, según cómo
termina el alineamiento en ella (algoritmo de Gotoh):

    M   emparejando a[i-1] con b[j-1]
    X   con a[i-1] frente a un hueco      (en la tabla se llega desde arriba)
    Y   con un hueco frente a b[j-1]      (se llega desde la izquierda)

Así se sabe si un hueco se abre (viene de M) o se alarga (viene de X o de Y). Con el
mismo coste para abrir y alargar, el resultado es el Needleman-Wunsch clásico.
"""

from dataclasses import dataclass, field
from typing import Mapping, Protocol

from .constants import BLOSUM62
from .validation import validate_dna, validate_protein

GAP = "-"
_NEG = -10 ** 9            # "menos infinito": una casilla imposible
_M, _X, _Y = 0, 1, 2


class ScoringScheme(Protocol):
    gap_open: int          # primera posición de un hueco
    gap_extend: int        # cada posición más del mismo hueco

    def pair(self, x: str, y: str) -> int: ...


@dataclass(frozen=True)
class Scoring:
    """Puntuación sencilla (ADN): todas las posiciones de un hueco cuestan lo mismo."""

    match: int = 1
    mismatch: int = -1
    gap: int = -2

    def pair(self, x: str, y: str) -> int:
        return self.match if x == y else self.mismatch

    @property
    def gap_open(self) -> int:
        return self.gap

    @property
    def gap_extend(self) -> int:
        return self.gap


@dataclass(frozen=True)
class MatrixScoring:
    """Puntuación con matriz de sustitución y huecos afines (proteínas).

    Un hueco de k posiciones suma gap_open + (k - 1) * gap_extend. Los valores por defecto
    son los de BLAST para proteínas: BLOSUM62, apertura 11 y extensión 1, es decir,
    -12 la primera posición y -1 cada una de las siguientes.
    """

    matrix: Mapping[str, Mapping[str, int]] = field(default_factory=lambda: BLOSUM62,
                                                     hash=False, repr=False)
    gap_open: int = -12
    gap_extend: int = -1

    def pair(self, x: str, y: str) -> int:
        return self.matrix[x][y]


# Tipo de secuencia -> (validación, puntuación por defecto)
KINDS = {
    "dna": (validate_dna, Scoring()),
    "protein": (validate_protein, MatrixScoring()),
}


@dataclass(frozen=True)
class Alignment:
    """Resultado de un alineamiento. Las dos cadenas tienen la misma longitud."""

    first: str
    second: str
    score: int
    kind: str = "dna"

    @property
    def length(self) -> int:
        return len(self.first)

    @property
    def matches(self) -> int:
        return sum(1 for a, b in zip(self.first, self.second) if a == b)

    @property
    def gaps(self) -> int:
        return sum(1 for a, b in zip(self.first, self.second) if GAP in (a, b))

    @property
    def mismatches(self) -> int:
        return self.length - self.matches - self.gaps

    @property
    def identity(self) -> float:
        """Porcentaje de posiciones idénticas sobre la longitud del alineamiento."""
        return self.matches / self.length * 100

    @property
    def markup(self) -> str:
        """Línea intermedia: '|' coincidencia, '.' sustitución, ' ' hueco."""
        return "".join(
            " " if GAP in (a, b) else "|" if a == b else "."
            for a, b in zip(self.first, self.second)
        )

    def pairs(self) -> list[tuple[int, int]]:
        """Qué posición de una secuencia queda frente a cuál de la otra (desde 0).

        Solo las columnas sin hueco. Sirve para relacionar dos proteínas aminoácido a
        aminoácido (p. ej. resaltar el mismo residuo en dos estructuras 3D):
            -ACG / TACG  ->  [(0, 1), (1, 2), (2, 3)]
        """
        result, i, j = [], 0, 0
        for a, b in zip(self.first, self.second):
            if a != GAP and b != GAP:
                result.append((i, j))
            i += a != GAP
            j += b != GAP
        return result


def align(first: str, second: str, scoring: ScoringScheme | None = None, *,
          kind: str = "dna") -> Alignment:
    """Alineamiento global óptimo de dos secuencias de ADN (kind="dna") o de proteínas
    (kind="protein"). Si no se indica puntuación, se usa la de ese tipo."""
    if kind not in KINDS:
        raise ValueError(f"Tipo de secuencia desconocido: {kind}")
    validate, default_scoring = KINDS[kind]
    a, b = validate(first), validate(second)
    aligned_a, aligned_b, score = _gotoh(a, b, scoring or default_scoring)
    return Alignment(aligned_a, aligned_b, score, kind)


def _best(m: int, x: int, y: int) -> tuple[int, int]:
    """La mayor de tres puntuaciones y de qué estado viene. En caso de empate prefiere
    M (emparejar), después X y después Y: así el resultado no depende del azar."""
    if m >= x and m >= y:
        return m, _M
    if x >= y:
        return x, _X
    return y, _Y


def _gotoh(a: str, b: str, scoring: ScoringScheme) -> tuple[str, str, int]:
    n, m = len(a), len(b)
    open_, extend = scoring.gap_open, scoring.gap_extend
    pair = scoring.pair

    # Solo se guardan dos filas de puntuaciones (la anterior y la actual), pero sí todas
    # las "flechas" para reconstruir el camino: 1 byte por casilla con los tres estados
    # empaquetados (2 bits cada uno). Para dos proteínas de 1.000 aminoácidos, ~1 MB.
    trace = [bytearray(m + 1) for _ in range(n + 1)]

    # Fila 0: alinear con nada = solo huecos (estado Y, desde la izquierda).
    m_prev = [_NEG] * (m + 1)
    x_prev = [_NEG] * (m + 1)
    y_prev = [_NEG] * (m + 1)
    m_prev[0] = 0
    for j in range(1, m + 1):
        y_prev[j] = open_ + (j - 1) * extend
        trace[0][j] = (_Y if j > 1 else _M) << 4

    for i in range(1, n + 1):
        m_cur = [_NEG] * (m + 1)
        x_cur = [_NEG] * (m + 1)
        y_cur = [_NEG] * (m + 1)
        # Columna 0: solo huecos (estado X, desde arriba).
        x_cur[0] = open_ + (i - 1) * extend
        trace[i][0] = (_X if i > 1 else _M) << 2
        row, ai = trace[i], a[i - 1]

        for j in range(1, m + 1):
            # M: emparejar a[i-1] con b[j-1], viniendo de la diagonal.
            best_m, from_m = _best(m_prev[j - 1], x_prev[j - 1], y_prev[j - 1])
            m_cur[j] = best_m + pair(ai, b[j - 1])
            # X: a[i-1] frente a un hueco. Desde X se alarga el hueco; desde M o Y se abre.
            best_x, from_x = _best(m_prev[j] + open_, x_prev[j] + extend, y_prev[j] + open_)
            x_cur[j] = best_x
            # Y: un hueco frente a b[j-1]. Desde Y se alarga; desde M o X se abre.
            best_y, from_y = _best(m_cur[j - 1] + open_, x_cur[j - 1] + open_,
                                   y_cur[j - 1] + extend)
            y_cur[j] = best_y
            row[j] = from_m | (from_x << 2) | (from_y << 4)

        m_prev, x_prev, y_prev = m_cur, x_cur, y_cur

    # Reconstrucción (traceback): desde la esquina final hacia atrás, siguiendo las flechas.
    score, state = _best(m_prev[m], x_prev[m], y_prev[m])
    aligned_a, aligned_b = [], []
    i, j = n, m
    while i > 0 or j > 0:
        arrows = trace[i][j]
        if state == _M:
            aligned_a.append(a[i - 1])
            aligned_b.append(b[j - 1])
            state = arrows & 3
            i, j = i - 1, j - 1
        elif state == _X:
            aligned_a.append(a[i - 1])
            aligned_b.append(GAP)
            state = (arrows >> 2) & 3
            i -= 1
        else:
            aligned_a.append(GAP)
            aligned_b.append(b[j - 1])
            state = (arrows >> 4) & 3
            j -= 1

    return "".join(reversed(aligned_a)), "".join(reversed(aligned_b)), score
