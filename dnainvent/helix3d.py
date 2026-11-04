"""Geometría 3D de la doble hélice de ADN (forma B, la más común en las células).

Inspirado en la idea de codegiovanni/DNA_helix (https://github.com/codegiovanni/DNA_helix):
describir cada hebra con una hélice paramétrica y girarla con una matriz de
rotación. El código es propio y los parámetros siguen las medidas reales del ADN.

La hélice paramétrica: el par de bases número i está a una altura i * RISE
sobre el eje Z y girado i * TWIST alrededor de él.

    x = RADIUS * cos(i * TWIST)
    y = RADIUS * sin(i * TWIST)
    z = i * RISE

Como el ángulo crece a la vez que la altura, la hélice es dextrógira (gira hacia
la derecha, como un tornillo normal), igual que el ADN B.

La segunda hebra no está en el lado opuesto (180°), sino a GROOVE grados.
Por un lado las hebras quedan más juntas (surco menor) y por el otro más
separadas (surco mayor): es donde muchas proteínas "leen" el ADN.

Las unidades son ángstroms (1 Å = 0,1 nanómetros).
"""

import math
from dataclasses import dataclass

from .constants import COMPLEMENT
from .errors import InvalidParameterError
from .validation import validate_dna

RISE = 3.38               # Å que sube la hélice por cada par de bases
BASE_PAIRS_PER_TURN = 10.5
TWIST = 360 / BASE_PAIRS_PER_TURN   # ≈ 34,3° de giro por par de bases
RADIUS = 10.0             # Å desde el eje hasta el esqueleto de azúcar-fosfato
GROOVE = 140.0            # ° entre las dos hebras, visto desde el eje

Point = tuple[float, float, float]


@dataclass(frozen=True)
class BasePair:
    """Un peldaño de la escalera: una base de cada hebra y dónde está su esqueleto.

    index      Posición en la secuencia (desde 0).
    base       Base de la hebra introducida (5'->3').
    partner    Base complementaria de la otra hebra.
    strand1    Punto del esqueleto de la hebra introducida.
    strand2    Punto del esqueleto de la hebra complementaria.
    """

    index: int
    base: str
    partner: str
    strand1: Point
    strand2: Point


def backbone_point(index: int, offset_degrees: float = 0.0) -> Point:
    """Punto del esqueleto del par 'index', girado 'offset_degrees' extra."""
    angle = math.radians(index * TWIST + offset_degrees)
    return (RADIUS * math.cos(angle), RADIUS * math.sin(angle), index * RISE)


def build_helix(dna: str, start: int = 0, count: int | None = None) -> list[BasePair]:
    """Pares de bases desde 'start' (como mucho 'count'), con sus coordenadas 3D.

    Las coordenadas dependen de la posición real en la secuencia, así que una
    ventana [start, start + count) encaja exactamente en la hélice completa.
    """
    dna = validate_dna(dna)
    if not 0 <= start < len(dna):
        raise InvalidParameterError("start_out_of_range", start=start, length=len(dna))
    end = len(dna) if count is None else min(len(dna), start + count)
    return [
        BasePair(i, dna[i], COMPLEMENT[dna[i]], backbone_point(i), backbone_point(i, GROOVE))
        for i in range(start, end)
    ]


def rotate_z(point: Point, theta: float, center: Point = (0.0, 0.0, 0.0)) -> Point:
    """Gira un punto 'theta' radianes alrededor de un eje paralelo a Z que pasa por 'center'.

    Es la matriz de rotación sobre Z (la misma idea que en DNA_helix):
        | cos θ  -sin θ  0 |
        | sin θ   cos θ  0 |
        |   0       0    1 |
    La web aplica esta misma fórmula 60 veces por segundo para que la hélice gire.
    """
    c, s = math.cos(theta), math.sin(theta)
    x, y, z = (p - q for p, q in zip(point, center))
    return (center[0] + c * x - s * y, center[1] + s * x + c * y, center[2] + z)
