"""
Nivel 2: un átomo dentro de la molécula (unidad: picómetro, 1 pm = 10⁻¹² m).

Se ve:
  - La nube de electrones del átomo, muestreada de |ψ|² con la física de AtomicWorld
    (atomicworld/scene/orbitals.py): orbitales hidrogenoides con la carga efectiva de
    Slater; en modo didáctico, los lóbulos p reales ocupados según la regla de Hund.
  - Las nubes de sus vecinos más cercanos: donde se solapan está el enlace covalente
    (electrones compartidos).
  - El resto de la molécula, de fondo.
  - El núcleo, a su tamaño real: ~0,003 pm. Desde fuera es invisible; hay que viajar por
    un "vacío" enorme hasta verlo aparecer.

Limitación: son orbitales de átomo aislado; en la molécula se deforman al enlazarse.
"""
import numpy as np

from atomicworld import atoms
from atomicworld.scene import orbitals
from dnainvent.scene.elements import atom_type
from dnainvent.scene.molecule import Molecule
from microscope.levels.base import Child, Focus, Level, Target
from microscope.levels.helix import atom_context, atom_title
from microscope.render.points import LineSet, PointCloud, PointStyle

SHELL_COLORS = {1: (1.00, 0.72, 0.30), 2: (0.30, 0.85, 1.00), 3: (0.78, 0.50, 1.00)}
REALISTIC_COLOR = np.array([0.55, 0.70, 1.00])
NUCLEUS_COLOR = (1.0, 0.45, 0.35, 1.0)
MAIN_POINTS = 90_000     # puntos totales de la nube principal (mín. por electrón: 6000)
NEIGHBOR_POINTS = 12_000
NEIGHBOR_RADIUS_NM = 0.32
EXIT_FRACTION = 0.8      # se sale a 0,8 × radio de van der Waals
ENTRY_FACTOR = 3.0       # se entra en el núcleo a 3 × su radio
_SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")


def config_string(shells) -> str:
    """'1s² 2s² 2p²' para el carbono."""
    return " ".join(f"{s.name}{str(s.electrons).translate(_SUP)}" for s in shells)


class AtomLevel(Level):
    depth = 1
    unit_m = 1e-12

    def __init__(self, molecule: Molecule, index: int, tr):
        super().__init__(tr)
        self.mol = molecule
        self.i = index
        self.type = atom_type(molecule.elements[index])
        self.vdw_pm = self.type.vdw_nm * 1000
        self.rel_pm = (molecule.positions - molecule.positions[index]) * 1000.0
        dist_nm = np.linalg.norm(self.rel_pm, axis=1) / 1000
        self.neighbors = [int(j) for j in np.argsort(dist_nm)
                          if 0 < dist_nm[j] < NEIGHBOR_RADIUS_NM]
        self.bonded = molecule.bonded_to(index)
        self.nucleus_r_pm = atoms.nuclear_radius_fm(self.type.isotope.a) / 1000
        self.clouds = {}
        self.lines = []

    @property
    def name(self):
        return atom_title(self.tr, self.mol, self.i)

    # ------------------------------------------------------------ construcción
    def _cloud(self, renderer, t, center, ppe, didactic, alpha, rng):
        pts, shells_n, _ = orbitals.sample_atom(t.number, t.shells, ppe, didactic, rng)
        if didactic:
            col = np.array([SHELL_COLORS[n] for n in shells_n])
        else:
            col = np.tile(REALISTIC_COLOR, (len(pts), 1))
        col = np.hstack([col, np.full((len(pts), 1), alpha)])
        seed = rng.random(len(pts)).astype("f4") * 1000
        return PointCloud(renderer, pts + center, col, np.full(len(pts), 0.9), seed)

    def build(self, renderer, didactic):
        rng = np.random.default_rng(self.i)
        t = self.type
        self.clouds["main"] = self._cloud(renderer, t, np.zeros(3),
                                          max(6000, MAIN_POINTS // t.electrons), didactic, 0.30, rng)
        if self.neighbors:
            parts = []
            for j in self.neighbors:
                tj = atom_type(self.mol.elements[j])
                parts.append(self._cloud(renderer, tj, self.rel_pm[j],
                                         max(800, NEIGHBOR_POINTS // tj.electrons), didactic, 0.22, rng))
            self.clouds["neighbors"] = parts

        # resto de la molécula como puntos de fondo
        others = [j for j in range(self.mol.n_atoms) if j != self.i and j not in self.neighbors]
        if others:
            types = [atom_type(self.mol.elements[j]) for j in others]
            col = np.array([[*(tj.color if didactic else REALISTIC_COLOR), 0.35] for tj in types])
            size = np.array([tj.vdw_nm * 1000 * 0.12 for tj in types])
            self.clouds["background"] = PointCloud(renderer, self.rel_pm[others], col, size,
                                                   np.array(others, dtype="f4"))

        # núcleo a escala real
        self.clouds["nucleus"] = PointCloud(renderer, np.zeros((1, 3)), np.array([NUCLEUS_COLOR]),
                                            np.array([self.nucleus_r_pm]), np.zeros(1))
        if didactic and self.bonded:
            a = np.zeros((len(self.bonded), 3))
            b = self.rel_pm[self.bonded]
            self.lines.append(LineSet(renderer, a, b, np.tile([1.0, 1.0, 1.0, 0.25], (len(a), 1))))

    def release(self):
        for v in self.clouds.values():
            for c in (v if isinstance(v, list) else [v]):
                c.release()
        for line in self.lines:
            line.release()
        self.clouds, self.lines = {}, []

    # ------------------------------------------------------------ dibujo
    def render(self, renderer, cam, didactic, fade):
        cloud_style = PointStyle(twinkle=0.5, min_px=1.5, max_px=5, softness=0.8)
        renderer.points(self.clouds["main"], cam, cloud_style, fade)
        for c in self.clouds.get("neighbors", []):
            renderer.points(c, cam, cloud_style, fade)
        if "background" in self.clouds:
            renderer.points(self.clouds["background"], cam,
                            PointStyle(min_px=1.5, max_px=300, softness=1.0), fade * 0.6)
        for line in self.lines:
            renderer.lines(line, cam, fade, 0.0)
        # El núcleo: en modo realista, tamaño real (invisible hasta acercarse mucho).
        # En modo didáctico, un marcador mínimo de 5 px para saber dónde está.
        renderer.points(self.clouds["nucleus"], cam,
                        PointStyle(min_px=5.0 if didactic else 0.5, max_px=900, softness=0.6), fade)

    # ------------------------------------------------------------ navegación
    def surface_distance(self, pos):
        d_nuc = np.linalg.norm(pos) - self.nucleus_r_pm
        d_nb = min((np.linalg.norm(pos - self.rel_pm[j]) for j in self.neighbors), default=np.inf)
        return float(max(min(d_nuc, d_nb), self.nucleus_r_pm * 0.5))

    def child_to_enter(self, pos, only=None):
        if only is not None and only != ("nucleus", self.i):
            return None
        if np.linalg.norm(pos) < ENTRY_FACTOR * self.nucleus_r_pm:
            from microscope.levels.nucleus import NucleusLevel
            mol, i = self.mol, self.i
            return Child(key=("nucleus", i), center=np.zeros(3), scale=1000.0,
                         create=lambda: NucleusLevel(self.type, self.tr,
                                                     context=lambda: atom_context(self.tr, mol, i)))
        return None

    def dive_target(self, pos, forward, didactic):
        tr, t = self.tr, self.type
        return Target(lambda: tr("microscope.autopilot.nucleus_target",
                                 of_element=tr.of(tr.element(t.number)), a=t.isotope.a),
                      lambda: np.zeros(3), key=("nucleus", self.i))

    def ascend_point(self, pos):
        r = np.linalg.norm(pos)
        direction = pos / r if r > 1e-9 else np.array([1.0, 0.0, 0.0])
        return direction * EXIT_FRACTION * self.vdw_pm * 1.5

    def should_exit(self, pos):
        r = np.linalg.norm(pos)
        if r > EXIT_FRACTION * self.vdw_pm:
            return True
        # si un vecino está más cerca que nuestro núcleo, salimos para entrar en él
        return any(np.linalg.norm(pos - self.rel_pm[j]) < r for j in self.bonded)

    def focus(self, pos, forward, didactic):
        tr, t = self.tr, self.type
        element = tr.element(t.number)
        r = float(np.linalg.norm(pos))
        to_nuc = -pos / max(r, 1e-12)
        if forward @ to_nuc > np.cos(np.radians(4)) or r < 5:
            title = tr("microscope.atom.nucleus_title", element=element, of_element=tr.of(element),
                       protons=tr.count(t.element.protons, "proton"),
                       neutrons=tr.count(t.isotope.neutrons, "neutron"))
            expl = []
            if didactic:
                ratio = self.vdw_pm / self.nucleus_r_pm
                expl = tr.lines("microscope.atom.nucleus_explanation",
                                r=tr.decimal(self.nucleus_r_pm * 1000, 1),
                                ratio=tr.number(int(round(ratio, -3))))
            return Focus(title, tr("microscope.atom.nucleus_distance", d=tr.decimal(r, 2)), expl,
                         distance=r)
        expl = []
        if didactic:
            zeff = ", ".join(f"{s.name} {tr.decimal(atoms.slater_zeff(t.element.protons, t.shells, s.n, s.l), 2)}"
                             for s in t.shells)
            expl = ([tr("microscope.atom.config", config=config_string(t.shells), zeff=zeff)]
                    + tr.lines("microscope.atom.cloud_explanation"))
        else:
            expl = tr.lines("microscope.atom.cloud_explanation_realistic")
        return Focus(tr("microscope.atom.cloud_title", electrons=tr.count(t.electrons, "electron")),
                     atom_context(tr, self.mol, self.i), expl,
                     pair_index=int(self.mol.pair_index[self.i]), distance=r)

    def time_info(self):
        return self.tr("microscope.atom.time")
