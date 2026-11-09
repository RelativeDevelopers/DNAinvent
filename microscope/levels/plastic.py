"""
Modo plástico: dónde va cada bola y cómo se mueve.

Es una ilustración, como un modelo de bolas de un aula:
  - Electrones (e⁻): cada uno da vueltas en un círculo del radio medio de su subnivel
    (el mismo que usan las nubes: radio de Slater), en un plano inclinado al azar. Es el
    modelo de Bohr (1913): útil para contar electrones y capas, pero los electrones no
    tienen trayectoria. Los textos del HUD lo dicen.
  - Protones (+) y neutrones (0): en las posiciones del núcleo de AtomicWorld, con un
    pequeño temblor (el movimiento de Fermi).
  - Quarks (u, d) y gluones (g): en las mismas posiciones que el modo didáctico.

Los electrones y los quarks son puntuales: sus bolas tienen un tamaño inventado y un
mínimo en píxeles para que siempre se vean.
"""
from dataclasses import dataclass

import numpy as np

from atomicworld.atoms import slater_zeff
from atomicworld.scene import orbitals

NUCLEON_BALL_FM = 0.75        # radio de la bola de un nucleón (radio de carga del protón:
                              # 0,84 fm; un poco menos para que se vean separados)
ELECTRON_BALL_FRACTION = 0.03  # radio de la bola del electrón / radio del subnivel más externo
ELECTRON_MIN_PX = 4.5
NUCLEON_MIN_PX = 2.0          # el núcleo visto desde el átomo: un puntito de colores
FERMI_WOBBLE_FM = 0.12        # temblor de los nucleones
OUTER_PERIOD_S = 7.0          # una vuelta del electrón más externo, en segundos de pantalla
INNER_PERIOD_S = 1.4          # los internos giran más deprisa, pero no tanto que parpadeen


@dataclass
class Orbits:
    radius: np.ndarray        # (N,) pm
    u: np.ndarray             # (N,3) los dos ejes del plano de cada órbita
    v: np.ndarray
    omega: np.ndarray         # (N,) rad/s
    phase: np.ndarray         # (N,)
    shell_n: np.ndarray       # (N,) capa de cada electrón

    @property
    def count(self) -> int:
        return len(self.radius)

    def positions(self, t: float, center=np.zeros(3)) -> np.ndarray:
        a = (self.omega * t + self.phase)[:, None]
        return center + self.radius[:, None] * (np.cos(a) * self.u + np.sin(a) * self.v)


def electron_orbits(z: int, shells, rng) -> Orbits:
    """Una órbita por electrón. Los electrones de un mismo subnivel se reparten en fase
    (y en planos distintos) para que no vayan en fila."""
    radius, u, v, omega, phase, shell_n = [], [], [], [], [], []
    radii = {(s.n, s.l): orbitals.mean_radius_pm(s.n, s.l, slater_zeff(z, shells, s.n, s.l))
             for s in shells}
    outer = max(radii.values())
    for sh in shells:
        r = radii[(sh.n, sh.l)]
        # 3.ª ley de Kepler (ω ∝ r^-1,5), entre los dos periodos de pantalla
        period = float(np.clip(OUTER_PERIOD_S * (r / outer) ** 1.5, INNER_PERIOD_S, OUTER_PERIOD_S))
        for e in range(sh.electrons):
            normal = rng.normal(size=3)
            normal /= np.linalg.norm(normal)
            a = np.cross(normal, [1.0, 0.0, 0.0] if abs(normal[0]) < 0.9 else [0.0, 1.0, 0.0])
            a /= np.linalg.norm(a)
            radius.append(r * rng.uniform(0.92, 1.08))
            u.append(a)
            v.append(np.cross(normal, a))
            omega.append(2 * np.pi / period * (1 if e % 2 == 0 else -1))   # espín ↑ / ↓
            phase.append(2 * np.pi * e / max(sh.electrons, 1) + rng.uniform(0, 0.6))
            shell_n.append(sh.n)
    return Orbits(np.array(radius), np.array(u), np.array(v), np.array(omega),
                  np.array(phase), np.array(shell_n))


def electron_ball_pm(z: int, shells) -> float:
    """Radio (inventado) de la bola de un electrón: proporcional al tamaño del átomo."""
    outer = max(orbitals.mean_radius_pm(s.n, s.l, slater_zeff(z, shells, s.n, s.l)) for s in shells)
    return ELECTRON_BALL_FRACTION * outer


def nucleon_kinds(kinds) -> list[str]:
    """['p', 'n', ...] -> ['proton', 'neutron', ...] (nombres de las imágenes)."""
    return ["proton" if k == "p" else "neutron" for k in kinds]


class Wobble:
    """Temblor de los nucleones alrededor de su sitio (tres senos por nucleón)."""

    def __init__(self, count: int, amplitude: float, rng):
        self.amplitude = amplitude
        self.freq = rng.uniform(2.0, 4.0, (count, 3))
        self.phase = rng.uniform(0, 2 * np.pi, (count, 3))

    def offsets(self, t: float) -> np.ndarray:
        return self.amplitude * np.sin(self.freq * t + self.phase)
