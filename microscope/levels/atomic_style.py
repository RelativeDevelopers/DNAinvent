"""
Modo atómico: los átomos de la molécula dibujados como el "Átomo en 3D" de AtomicWorld.

Qué lo hace distinto del modo didáctico:
  - Colores por capa (naranja n = 1, cian n = 2, violeta n = 3), no por elemento.
  - Cada subnivel recibe más puntos cuanto mayor es (orbitals.didactic_points): con el
    mismo número de puntos por electrón, los electrones internos (pocos, apretados) brillan
    muchísimo y la capa externa (la que da tamaño y química al átomo) casi no se ve.
  - Brillo y tamaño de punto por subnivel (orbitals.didactic_style): las capas internas
    más tenues, las externas con puntos algo mayores para que formen una nube.

Toda la física sale de AtomicWorld (atomicworld/scene/orbitals.py): con la misma semilla,
el átomo de oxígeno se ve igual en las dos apps.
"""
import numpy as np

from atomicworld.scene import orbitals

SHELL_COLORS = {1: (1.00, 0.72, 0.30), 2: (0.30, 0.85, 1.00), 3: (0.78, 0.50, 1.00),
                4: (0.40, 1.00, 0.55), 5: (1.00, 0.45, 0.70), 6: (1.00, 0.95, 0.40),
                7: (0.50, 0.60, 1.00)}
# un color por orbital en la vista de un subnivel (tecla O): se distingue qué lóbulo es de cuál
ORBITAL_COLORS = [(1.00, 0.55, 0.25), (0.30, 0.85, 1.00), (0.55, 1.00, 0.40), (1.00, 0.40, 0.75),
                  (1.00, 0.92, 0.35), (0.65, 0.55, 1.00), (0.35, 1.00, 0.85)]
AXIS_COLORS = {"x": (1.0, 0.45, 0.45), "y": (0.45, 1.0, 0.45), "z": (0.45, 0.65, 1.0)}
ATOM_BUDGET = 220_000          # puntos del átomo principal (los mismos que en AtomicWorld)
MAX_ALPHA = 0.38               # brillo del subnivel más externo en orbitals.didactic_style
HELIX_ALPHA = 0.30             # brillo de la capa externa en la hélice
NUCLEUS_COLOR = (1.0, 0.45, 0.35)


def shell_colors(shells_n: np.ndarray) -> np.ndarray:
    palette = np.array([SHELL_COLORS[n] for n in sorted(SHELL_COLORS)])
    return palette[shells_n - 1]


def atom_cloud(atom_type, budget: int, rng):
    """Nube de un átomo como en AtomicWorld.

    Devuelve (posiciones_pm (N,3), rgb (N,3), alpha (N,), tamaño_pm (N,))."""
    z, shells = atom_type.element.protons, atom_type.shells
    ppe = orbitals.didactic_points(z, shells, budget)
    pts, shells_n, shells_l = orbitals.sample_atom(z, shells, ppe, True, rng)
    alpha = np.empty(len(pts), dtype="f4")
    size = np.empty(len(pts), dtype="f4")
    for (n, l), (a, sz) in orbitals.didactic_style(z, shells).items():
        mask = (shells_n == n) & (shells_l == l)
        alpha[mask], size[mask] = a, sz
    return pts, shell_colors(shells_n), alpha, size


def molecule_cloud(mol, atom_type, total_points: int, rng):
    """Todos los átomos de una molécula (hélice o proteína), en nanómetros.

    Cada átomo recibe el mismo número de puntos (como en AtomicWorld, donde el hidrógeno y
    el uranio se ven igual de llenos), repartidos entre sus subniveles como en atom_cloud.
    Devuelve (pos, color RGBA, tamaño, semilla), listos para un PointCloud."""
    per_atom = int(np.clip(total_points / max(mol.n_atoms, 1), 60, 1500))
    templates = {}
    for sym in set(mol.elements):
        t = atom_type(sym)
        pts, rgb, alpha, size_pm = atom_cloud(t, per_atom * 40, rng)
        # mismo reparto de brillo que en el átomo, escalado al de la hélice; puntos mucho
        # mayores que en el átomo: aquí cada uno representa muchos
        a = HELIX_ALPHA * alpha / MAX_ALPHA
        size_nm = 0.008 + 0.0012 * size_pm
        templates[sym] = (pts / 1000.0, rgb, a, size_nm)

    pos, col, size, seed = [], [], [], []
    for i, sym in enumerate(mol.elements):
        src, rgb, a, sz = templates[sym]
        k = rng.integers(0, len(src), per_atom)
        pos.append(mol.positions[i] + src[k])
        col.append(np.hstack([rgb[k], a[k, None]]))
        size.append(sz[k])
        seed.append(np.full(per_atom, float(i)))
    return np.vstack(pos), np.vstack(col), np.concatenate(size), np.concatenate(seed)
