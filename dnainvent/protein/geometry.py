"""Geometría 3D mínima en Python puro (vectores como tuplas), sin NumPy.

La pieza clave es place(): colocar un átomo a partir de sus "coordenadas internas".
En vez de decir "el hidrógeno está en (x, y, z)", se dice:

    está a r Å del átomo c,
    formando un ángulo θ con b-c,
    y girado un ángulo φ (diedro) alrededor del eje b-c respecto a a.

        a
         \\
          b ---- c
                  \\  r, θ
                   H        (φ: cuánto gira H alrededor de b-c, visto desde a)

Estas tres medidas no cambian aunque la proteína se mueva o gire, así que sirven para
guardar "dónde va cada hidrógeno" una sola vez por aminoácido.
"""

import math

Vec = tuple[float, float, float]


def sub(a: Vec, b: Vec) -> Vec:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def add(a: Vec, b: Vec) -> Vec:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def scale(a: Vec, k: float) -> Vec:
    return (a[0] * k, a[1] * k, a[2] * k)


def dot(a: Vec, b: Vec) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a: Vec, b: Vec) -> Vec:
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def norm(a: Vec) -> float:
    return math.sqrt(dot(a, a))


def unit(a: Vec) -> Vec:
    return scale(a, 1.0 / norm(a))


def distance(a: Vec, b: Vec) -> float:
    return norm(sub(a, b))


def angle(a: Vec, b: Vec, c: Vec) -> float:
    """Ángulo a-b-c en grados (vértice en b)."""
    u, v = unit(sub(a, b)), unit(sub(c, b))
    return math.degrees(math.acos(max(-1.0, min(1.0, dot(u, v)))))


def dihedral(a: Vec, b: Vec, c: Vec, d: Vec) -> float:
    """Ángulo diedro a-b-c-d en grados (-180, 180]: el giro de d respecto a a alrededor de b-c."""
    b0, b1, b2 = sub(a, b), sub(c, b), sub(d, c)
    b1 = unit(b1)
    v = sub(b0, scale(b1, dot(b0, b1)))
    w = sub(b2, scale(b1, dot(b2, b1)))
    return math.degrees(math.atan2(dot(cross(b1, v), w), dot(v, w)))


def place(a: Vec, b: Vec, c: Vec, r: float, theta: float, phi: float) -> Vec:
    """Posición de d tal que |c-d| = r, ángulo b-c-d = theta y diedro a-b-c-d = phi.

    Método NeRF (Natural Extension Reference Frame, Parsons et al. 2005): se construye un
    sistema de ejes en c a partir de a, b y c, y se coloca d en él.
    """
    theta, phi = math.radians(theta), math.radians(phi)
    bc = unit(sub(c, b))
    n = unit(cross(sub(b, a), bc))
    m = cross(n, bc)
    local = (-r * math.cos(theta), r * math.sin(theta) * math.cos(phi),
             r * math.sin(theta) * math.sin(phi))
    return add(c, add(add(scale(bc, local[0]), scale(m, local[1])), scale(n, local[2])))
