"""Superponer dos estructuras: el giro y desplazamiento que mejor encaja una sobre otra.

Cada estructura llega con sus coordenadas en una posición y orientación arbitrarias. Para
compararlas hay que moverlas hasta que encajen lo mejor posible: el movimiento rígido
(rotación R + traslación t) que hace mínima la distancia media entre átomos equivalentes.

Lo que queda de distancia después de superponer es el RMSD (raíz de la media de los
cuadrados de las distancias), en Å:

    < 1 Å   prácticamente la misma forma
    1-2 Å   la misma forma con pequeñas diferencias
    > 3 Å   formas claramente distintas (o mal superpuestas)

Método de Horn (1987, cuaterniones): la mejor rotación es el vector propio del mayor valor
propio de una matriz 4x4 que se construye con las coordenadas. Da lo mismo que el método de
Kabsch (1976), pero sin descomposición SVD, así que cabe en Python puro: el vector propio
se calcula con el método de Jacobi (giros sucesivos que anulan los elementos fuera de la
diagonal), que para una matriz 4x4 converge en unas pocas vueltas.
"""

import math
from dataclasses import dataclass

from .geometry import Vec, add, scale, sub

Matrix = tuple[Vec, Vec, Vec]


@dataclass(frozen=True)
class Superposition:
    rotation: Matrix          # R (por filas)
    translation: Vec          # t
    rmsd: float               # Å, después de superponer
    count: int                # cuántos puntos se usaron

    def apply(self, point: Vec) -> Vec:
        """R · p + t: lleva un punto de la estructura móvil a la de referencia."""
        r = self.rotation
        return (sum(r[0][k] * point[k] for k in range(3)) + self.translation[0],
                sum(r[1][k] * point[k] for k in range(3)) + self.translation[1],
                sum(r[2][k] * point[k] for k in range(3)) + self.translation[2])


def superpose(mobile: list[Vec], reference: list[Vec]) -> Superposition:
    """Movimiento rígido que lleva 'mobile' encima de 'reference' (mismo número de puntos,
    en el mismo orden: el punto i de una corresponde al i de la otra)."""
    if len(mobile) != len(reference) or len(mobile) < 3:
        raise ValueError("Hacen falta al menos 3 parejas de puntos")
    n = len(mobile)
    c_mob, c_ref = _centroid(mobile), _centroid(reference)
    p = [sub(x, c_mob) for x in mobile]
    q = [sub(x, c_ref) for x in reference]

    # Matriz de correlación S[a][b] = Σ p_a · q_b
    s = [[sum(pi[a] * qi[b] for pi, qi in zip(p, q)) for b in range(3)] for a in range(3)]
    (sxx, sxy, sxz), (syx, syy, syz), (szx, szy, szz) = s
    horn = [
        [sxx + syy + szz, syz - szy, szx - sxz, sxy - syx],
        [syz - szy, sxx - syy - szz, sxy + syx, szx + sxz],
        [szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy],
        [sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz],
    ]
    values, vectors = _jacobi(horn)
    best = max(range(4), key=lambda k: values[k])
    rotation = _rotation_from_quaternion([vectors[i][best] for i in range(4)])

    moved_centroid = _rotate(rotation, c_mob)
    translation = sub(c_ref, moved_centroid)
    result = Superposition(rotation, translation, 0.0, n)
    squared = sum(_dist2(result.apply(a), b) for a, b in zip(mobile, reference))
    return Superposition(rotation, translation, math.sqrt(squared / n), n)


def rmsd(a: list[Vec], b: list[Vec]) -> float:
    """RMSD sin mover nada (las coordenadas ya están superpuestas)."""
    return math.sqrt(sum(_dist2(x, y) for x, y in zip(a, b)) / len(a))


def _centroid(points: list[Vec]) -> Vec:
    total = (0.0, 0.0, 0.0)
    for point in points:
        total = add(total, point)
    return scale(total, 1.0 / len(points))


def _dist2(a: Vec, b: Vec) -> float:
    d = sub(a, b)
    return d[0] * d[0] + d[1] * d[1] + d[2] * d[2]


def _rotate(r: Matrix, v: Vec) -> Vec:
    return tuple(sum(r[i][k] * v[k] for k in range(3)) for i in range(3))


def _rotation_from_quaternion(q: list[float]) -> Matrix:
    w, x, y, z = (c / math.sqrt(sum(c * c for c in q)) for c in q)
    return (
        (w * w + x * x - y * y - z * z, 2 * (x * y - w * z), 2 * (x * z + w * y)),
        (2 * (x * y + w * z), w * w - x * x + y * y - z * z, 2 * (y * z - w * x)),
        (2 * (x * z - w * y), 2 * (y * z + w * x), w * w - x * x - y * y + z * z),
    )


def _jacobi(matrix: list[list[float]], sweeps: int = 50) -> tuple[list[float], list[list[float]]]:
    """Valores y vectores propios (en columnas) de una matriz simétrica pequeña."""
    n = len(matrix)
    a = [row[:] for row in matrix]
    v = [[float(i == j) for j in range(n)] for i in range(n)]
    for _ in range(sweeps):
        off = sum(a[i][j] ** 2 for i in range(n) for j in range(n) if i != j)
        if off < 1e-22 * (1 + sum(a[i][i] ** 2 for i in range(n))):
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                if abs(a[p][q]) < 1e-300:
                    continue
                theta = (a[q][q] - a[p][p]) / (2 * a[p][q])
                t = math.copysign(1.0, theta) / (abs(theta) + math.sqrt(theta * theta + 1))
                c = 1 / math.sqrt(t * t + 1)
                s = t * c
                for k in range(n):                       # columnas p y q
                    akp, akq = a[k][p], a[k][q]
                    a[k][p], a[k][q] = c * akp - s * akq, s * akp + c * akq
                for k in range(n):                       # filas p y q
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k], a[q][k] = c * apk - s * aqk, s * apk + c * aqk
                for k in range(n):
                    vkp, vkq = v[k][p], v[k][q]
                    v[k][p], v[k][q] = c * vkp - s * vkq, s * vkp + c * vkq
    return [a[i][i] for i in range(n)], v
