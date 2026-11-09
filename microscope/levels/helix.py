"""
Nivel 1: la doble hélice completa (unidad: nanómetro).

Cada átomo se dibuja como su nube de electrones (los mismos orbitales que se ven en el
nivel atómico, con menos puntos). La forma de la molécula emerge sola de la suma de las
nubes, como en la realidad.

En el modo atómico, cada átomo se dibuja como en AtomicWorld (levels/atomic_style.py):
colores por capa y todas las capas visibles. En el modo plástico, la hélice se ve como en
el didáctico: las bolas de plástico aparecen al entrar en un átomo.
"""
import numpy as np

from atomicworld.scene import orbitals
from dnainvent.scene.elements import DNA_ELEMENTS, atom_type
from dnainvent.scene.molecule import Molecule
from microscope.levels.atomic_style import NUCLEUS_COLOR, molecule_cloud
from microscope.levels.base import Child, Focus, Level, Target, aimed
from microscope.render.points import LineSet, PointCloud, PointStyle

REALISTIC_COLOR = np.array([0.55, 0.70, 1.00])
TARGET_POINTS = 350_000
ENTER_FRACTION = 0.45          # entras en un átomo a 0,45 × su radio de van der Waals
VIB_AMP_NM = 0.006             # vibración térmica (~0,05 Å), ver time_info()
HBOND_COLOR = (1.0, 0.9, 0.3, 0.8)


def atom_context(tr, mol: Molecule, i: int) -> str:
    """'Adenina (A) · par 3 · hebra 1' o 'Grupo fosfato del nucleótido G · par 5 · hebra 2'.
    En una proteína: 'esqueleto · Glutamato (Glu) 7' (ver levels/protein.py)."""
    if getattr(mol, "kind", "dna") == "protein":
        from microscope.levels.protein import protein_atom_context
        return protein_atom_context(tr, mol, i)
    part, base = mol.part[i], mol.resn[i]
    if part == "base":
        what = tr("microscope.molecule.base", name=tr(f"microscope.bases.{base}"), base=base)
    else:
        what = tr("microscope.molecule.part_of", part=tr(f"microscope.parts.{part}"), base=base)
    return tr("microscope.molecule.where", what=what, pair=int(mol.pair_index[i]) + 1,
              strand=int(mol.strand[i]))


def atom_title(tr, mol: Molecule, i: int) -> str:
    """'Átomo de nitrógeno (N7)', 'Nitrogen atom (N7)', "Àtom de nitrogen (N7)"."""
    element = tr.element(atom_type(mol.elements[i]).number)
    return tr("microscope.molecule.atom", element=element, of_element=tr.of(element),
              name=mol.names[i])


class HelixLevel(Level):
    depth = 0
    unit_m = 1e-9

    def __init__(self, molecule: Molecule, tr):
        super().__init__(tr)
        self.mol = molecule
        self.vdw = molecule.radius_vdw()
        self.clouds = []
        self.lines = []

    @property
    def name(self):
        return self.tr("microscope.levels.helix")

    # ------------------------------------------------------------ construcción
    def build(self, renderer, didactic: bool):
        if self.mode == "atomic":
            return self._build_atomic(renderer)
        mol = self.mol
        rng = np.random.default_rng(1)
        types = {s: atom_type(s) for s in DNA_ELEMENTS}
        total_e = sum(types[e].electrons for e in mol.elements)
        ppe = int(np.clip(TARGET_POINTS / total_e, 4, 60))

        # una nube "de muestra" por elemento (pm -> nm); cada átomo toma puntos al azar de ella
        cache = {s: orbitals.sample_atom(t.number, t.shells, ppe * 40, didactic, rng)[0] / 1000.0
                 for s, t in types.items()}

        pos, col, size, seed = [], [], [], []
        for i, sym in enumerate(mol.elements):
            t, src = types[sym], cache[sym]
            n = t.electrons * ppe
            pos.append(mol.positions[i] + src[rng.integers(0, len(src), n)])
            c = np.array(t.color) if didactic else REALISTIC_COLOR
            col.append(np.tile([*c, 0.16], (n, 1)))
            size.append(np.full(n, 0.010))
            seed.append(np.full(n, float(i)))
        self.clouds.append(PointCloud(renderer, np.vstack(pos), np.vstack(col),
                                      np.concatenate(size), np.concatenate(seed)))

        if didactic:
            # "Núcleos" marcados como esferas pequeñas del color del elemento (no es a escala:
            # el núcleo real es 100.000 veces más pequeño que el átomo)
            cpk = np.array([[*types[e].color, 0.9] for e in mol.elements])
            self.clouds.append(PointCloud(renderer, mol.positions, cpk,
                                          self.vdw * 0.22, np.arange(mol.n_atoms, dtype="f4")))
            a, b = mol.positions[mol.bonds[:, 0]], mol.positions[mol.bonds[:, 1]]
            ca = np.array([[*types[mol.elements[i]].color, 0.55] for i in mol.bonds[:, 0]])
            cb = np.array([[*types[mol.elements[i]].color, 0.55] for i in mol.bonds[:, 1]])
            self.lines.append((LineSet(renderer, a, b, ca, cb), 0.0))
            if len(mol.hbonds):
                ha, hb = mol.positions[mol.hbonds[:, 0]], mol.positions[mol.hbonds[:, 1]]
                self.lines.append((LineSet(renderer, ha, hb, np.tile(HBOND_COLOR, (len(ha), 1))), 0.03))

    def _build_atomic(self, renderer):
        """Modo atómico: cada átomo como en AtomicWorld; solo los puentes de hidrógeno."""
        mol = self.mol
        self.clouds.append(PointCloud(renderer, *molecule_cloud(mol, atom_type, TARGET_POINTS,
                                                                np.random.default_rng(1))))
        # los núcleos, como en AtomicWorld: un puntito rojo (no a escala) para ver la estructura
        col = np.tile([*NUCLEUS_COLOR, 0.7], (mol.n_atoms, 1))
        self.clouds.append(PointCloud(renderer, mol.positions, col, np.full(mol.n_atoms, 0.012),
                                      np.arange(mol.n_atoms, dtype="f4")))
        if len(mol.hbonds):
            ha, hb = mol.positions[mol.hbonds[:, 0]], mol.positions[mol.hbonds[:, 1]]
            self.lines.append((LineSet(renderer, ha, hb, np.tile(HBOND_COLOR, (len(ha), 1))), 0.03))

    def release(self):
        for c in self.clouds:
            c.release()
        for line, _ in self.lines:
            line.release()
        self.clouds, self.lines = [], []

    # ------------------------------------------------------------ dibujo
    def render(self, renderer, cam, didactic, fade):
        style = PointStyle(vib_amp=VIB_AMP_NM, min_px=1.2, max_px=400, softness=1.0)
        renderer.points(self.clouds[0], cam, style, fade)
        if self.mode == "atomic":
            renderer.points(self.clouds[1], cam, PointStyle(vib_amp=VIB_AMP_NM, min_px=1.5,
                                                            max_px=8, softness=0.6), fade)
        elif didactic and len(self.clouds) > 1:
            renderer.points(self.clouds[1], cam,
                            PointStyle(vib_amp=VIB_AMP_NM, min_px=2.0, max_px=200,
                                       softness=0.35), fade * 0.9)
        for lines, dash in self.lines:
            renderer.lines(lines, cam, fade, dash)

    # ------------------------------------------------------------ navegación
    def _dist(self, pos):
        return np.linalg.norm(self.mol.positions - pos, axis=1)

    def surface_distance(self, pos):
        return float(self._dist(pos).min())

    def child_to_enter(self, pos, only=None):
        d = self._dist(pos)
        i = only[1] if only is not None else int(np.argmin(d))
        if d[i] < ENTER_FRACTION * self.vdw[i]:
            from microscope.levels.atom import AtomLevel
            return Child(key=("atom", i), center=self.mol.positions[i].astype(float),
                         scale=1000.0, create=lambda: AtomLevel(self.mol, i, self.tr))
        return None

    def dive_target(self, pos, forward, didactic):
        mol = self.mol
        # 1) lo que hay en la mira (aunque sea un hidrógeno: lo ha elegido el usuario)
        i = aimed(mol.positions, pos, forward, 5)
        if i is None:
            # 2) si no, el átomo pesado (no H) más cercano y razonablemente de frente
            v = mol.positions - pos
            d = np.linalg.norm(v, axis=1)
            cos = (v @ forward) / np.maximum(d, 1e-12)
            heavy = np.array([e != "H" for e in mol.elements])
            i = int(np.argmin(np.where(heavy, d * (2.0 - cos), np.inf)))
        p = mol.positions[i].astype(float)
        tr = self.tr
        return Target(lambda: tr("microscope.autopilot.atom_target", atom=atom_title(tr, mol, i),
                                 pair=int(mol.pair_index[i]) + 1),
                      lambda: p, key=("atom", i))

    def focus(self, pos, forward, didactic):
        mol, tr = self.mol, self.tr
        i = aimed(mol.positions, pos, forward, 5)
        if i is None:
            if self._dist(pos).min() > 6:
                return self._overview(didactic)
            return None
        return self.atom_focus(i, float(self._dist(pos)[i]), didactic)

    def _overview(self, didactic):
        mol, tr = self.mol, self.tr
        shown = mol.sequence[:24] + ("…" if len(mol.sequence) > 24 else "")
        if mol.start == 0 and len(mol.sequence) == mol.length:
            detail = tr("microscope.helix.pairs", pairs=len(mol.sequence), atoms=mol.n_atoms)
        else:
            detail = tr("microscope.helix.window", first=mol.start + 1, last=mol.end,
                        length=mol.length, atoms=mol.n_atoms)
        expl = tr.lines("microscope.helix.explanation") if didactic else []
        if self.mode in ("atomic", "plastic"):
            expl = expl + tr.lines(f"microscope.helix.explanation_{self.mode}")
        return Focus(tr("microscope.helix.title", sequence=shown), detail, expl)

    def atom_focus(self, i, dist, didactic):
        mol, tr = self.mol, self.tr
        expl = tr.lines(f"microscope.helix.explanation_{mol.part[i]}") if didactic else []
        return Focus(atom_title(tr, mol, i), atom_context(tr, mol, i), expl,
                     pair_index=int(mol.pair_index[i]), distance=dist)

    def time_info(self):
        return self.tr("microscope.helix.time")
