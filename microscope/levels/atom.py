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

Modos propios de este nivel:
  - atómico: el átomo como el "Átomo en 3D" de AtomicWorld (levels/atomic_style.py), con
    220.000 puntos repartidos para que se vean todas las capas, y las teclas O (un subnivel
    cada vez) e I (un orbital p cada vez) con los ejes x, y, z.
  - plástico: los electrones como bolas verdes (e⁻) dando vueltas, y el núcleo de bolas
    rojas (+) y azules (0) a escala real (levels/plastic.py).
"""
import numpy as np

from atomicworld import atoms
from atomicworld.scene import orbitals
from dnainvent.scene.elements import atom_type
from dnainvent.scene.molecule import Molecule
from atomicworld.scene import nucleus as nuc
from microscope.levels import atomic_style, plastic
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
ATOMIC_EXIT_FRACTION = 2.2   # modo atómico: se puede retroceder hasta ver el átomo entero
ENTRY_FACTOR = 3.0       # se entra en el núcleo a 3 × su radio
ATOMIC_NEIGHBOR_BUDGET = 30_000     # modo atómico: puntos de cada vecino
DETAIL_POINTS = 150_000             # modo atómico: puntos del subnivel elegido con O / I
PLASTIC_CLOUD_FADE = 0.25           # modo plástico: la nube de fondo, tenue
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
        self.shells = self.type.shells
        self.view = None            # modo atómico: None = todos; i = solo el subnivel shells[i]
        self.orbital = None         # con un subnivel p elegido: None = todos; k = solo uno
        self.renderer = None
        self.electrons = []         # modo plástico: [(Orbits, centro)] del átomo y sus vecinos
        self.balls = None           # modo plástico: SpriteSet con electrones y nucleones
        self.time = 0.0

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
        self.renderer = renderer
        if self.mode == "atomic":
            return self._build_atomic(renderer)
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
        self._background(renderer, didactic)

        # núcleo a escala real
        self.clouds["nucleus"] = PointCloud(renderer, np.zeros((1, 3)), np.array([NUCLEUS_COLOR]),
                                            np.array([self.nucleus_r_pm]), np.zeros(1))
        if didactic and self.bonded:
            self._bond_lines(renderer)
        if self.mode == "plastic":
            self._build_plastic(renderer)

    def _bond_lines(self, renderer):
        a = np.zeros((len(self.bonded), 3))
        b = self.rel_pm[self.bonded]
        self.lines.append(LineSet(renderer, a, b, np.tile([1.0, 1.0, 1.0, 0.25], (len(a), 1))))

    def _background(self, renderer, didactic):
        others = [j for j in range(self.mol.n_atoms) if j != self.i and j not in self.neighbors]
        if others:
            types = [atom_type(self.mol.elements[j]) for j in others]
            col = np.array([[*(tj.color if didactic else REALISTIC_COLOR), 0.35] for tj in types])
            size = np.array([tj.vdw_nm * 1000 * 0.12 for tj in types])
            self.clouds["background"] = PointCloud(renderer, self.rel_pm[others], col, size,
                                                   np.array(others, dtype="f4"))

    # ------------------------------------------------------------ modo atómico
    def _atomic_cloud(self, renderer, t, center, budget, alpha_scale, rng):
        pts, rgb, alpha, size = atomic_style.atom_cloud(t, budget, rng)
        col = np.hstack([rgb, (alpha * alpha_scale)[:, None]])
        seed = rng.random(len(pts)).astype("f4") * 1000
        return PointCloud(renderer, pts + center, col, size, seed)

    def _build_atomic(self, renderer):
        """El átomo como en AtomicWorld; los vecinos igual pero más tenues."""
        rng = np.random.default_rng(self.i)
        self.clouds["main"] = self._atomic_cloud(renderer, self.type, np.zeros(3),
                                                 atomic_style.ATOM_BUDGET, 1.0, rng)
        if self.neighbors:
            self.clouds["neighbors"] = [
                self._atomic_cloud(renderer, atom_type(self.mol.elements[j]), self.rel_pm[j],
                                   ATOMIC_NEIGHBOR_BUDGET, 0.45, rng) for j in self.neighbors]
        self._background(renderer, True)
        self.clouds["nucleus"] = PointCloud(renderer, np.zeros((1, 3)), np.array([NUCLEUS_COLOR]),
                                            np.array([self.nucleus_r_pm]), np.zeros(1))
        if self.bonded:
            self._bond_lines(renderer)
        self._build_detail()

    def cycle_view(self):
        """Tecla O. Todos → el subnivel más externo → ... → el 1s → todos otra vez."""
        order = [None] + list(range(len(self.shells) - 1, -1, -1))
        self.view = order[(order.index(self.view) + 1) % len(order)]
        self.orbital = None
        self._build_detail()
        return self.view

    def cycle_orbital(self):
        """Tecla I, con un subnivel p elegido: sus orbitales ocupados de uno en uno."""
        if self.view is None or self.shells[self.view].l == 0:
            return False
        sh = self.shells[self.view]
        occupied = [k for k, e in enumerate(orbitals.hund_occupancy(sh.l, sh.electrons)) if e]
        order = [None] + occupied
        self.orbital = order[(order.index(self.orbital) + 1) % len(order)]
        self._build_detail()
        return True

    def subshell_label(self, i) -> str:
        sh = self.shells[i]
        return f"{sh.name}{str(sh.electrons).translate(_SUP)}"

    def orbital_label(self) -> str:
        sh = self.shells[self.view]
        return f"{sh.n}{orbitals.orbital_names(sh.l)[self.orbital]}"

    def _zeff(self, sh) -> float:
        return atoms.slater_zeff(self.type.element.protons, self.shells, sh.n, sh.l)

    def axis_length(self) -> float:
        sh = self.shells[self.view]
        return 2.0 * orbitals.mean_radius_pm(sh.n, sh.l, self._zeff(sh))

    def axis_labels(self):
        """Puntas de los ejes (para escribir x, y, z en el HUD), o [] si no se muestran."""
        if "axes" not in self.clouds:
            return []
        half = self.axis_length() * 1.05
        return [(np.array(v, dtype=float) * half, k) for k, v in
                (("x", (1, 0, 0)), ("y", (0, 1, 0)), ("z", (0, 0, 1)))]

    def _build_detail(self):
        """El subnivel elegido con O, generado aparte con muchos puntos (en la vista de
        todos los electrones los internos llevan pocos y solos casi no se verían)."""
        for key in ("detail", "axes"):
            old = self.clouds.pop(key, None)
            if old:
                old.release()
        if self.view is None or self.renderer is None or self.mode != "atomic":
            return
        rng = np.random.default_rng(1000 + self.view)
        sh = self.shells[self.view]
        shown = sh.electrons if self.orbital is None else \
            orbitals.hund_occupancy(sh.l, sh.electrons)[self.orbital]
        pts, owner = orbitals.sample_subshell(self.type.element.protons, self.shells, self.view,
                                              DETAIL_POINTS // shown, True, rng, only=self.orbital)
        if sh.l == 0:
            rgb = np.tile(atomic_style.SHELL_COLORS[sh.n], (len(pts), 1))
        else:
            rgb = np.array(atomic_style.ORBITAL_COLORS)[owner % len(atomic_style.ORBITAL_COLORS)]
        col = np.hstack([rgb, np.full((len(pts), 1), 0.4)])
        self.clouds["detail"] = PointCloud(self.renderer, pts, col, np.full(len(pts), 0.9),
                                           rng.random(len(pts)).astype("f4") * 1000)
        half = self.axis_length()
        a = np.array([[-half, 0, 0], [0, -half, 0], [0, 0, -half]])
        colors = np.array([[*atomic_style.AXIS_COLORS[k], 0.55] for k in "xyz"])
        self.clouds["axes"] = LineSet(self.renderer, a, -a, colors)

    # ------------------------------------------------------------ modo plástico
    def _build_plastic(self, renderer):
        """Electrones (e⁻) del átomo y de sus vecinos, y sus núcleos con protones y neutrones."""
        from microscope.render.sprites import SpriteSet
        rng = np.random.default_rng(500 + self.i)
        self.electrons, radius, kinds, nucleons = [], [], [], []
        for j, center in [(self.i, np.zeros(3))] + [(j, self.rel_pm[j]) for j in self.neighbors]:
            t = atom_type(self.mol.elements[j])
            orbits = plastic.electron_orbits(t.element.protons, t.shells, rng)
            self.electrons.append((orbits, center))
            radius += [plastic.electron_ball_pm(t.element.protons, t.shells)] * orbits.count
            kinds += ["electron"] * orbits.count
            # núcleo a escala real (fm -> pm), con la misma semilla que el nivel del núcleo
            A = t.isotope.a
            nrng = np.random.default_rng(t.element.number * 1000 + A)
            centers = nuc.pack_nucleons(A, atoms.nuclear_radius_fm(A), nrng)
            nucleons.append(center + centers / 1000.0)
            radius += [plastic.NUCLEON_BALL_FM / 1000.0] * A
            kinds += plastic.nucleon_kinds(nuc.assign_kinds(t.element.protons, t.isotope.neutrons, nrng))
        self.nucleon_pos = np.vstack(nucleons)
        # electrones con un mínimo de 7 px; el núcleo, un puntito de colores que se abre en
        # bolas al acercarse
        min_px = np.where(np.array(kinds) == "electron", plastic.ELECTRON_MIN_PX,
                          plastic.NUCLEON_MIN_PX)
        self.balls = SpriteSet(renderer, self._ball_positions(self.time), np.array(radius), kinds,
                               min_px=min_px)

    def _ball_positions(self, t):
        electrons = [o.positions(t, c) for o, c in self.electrons]
        return np.vstack(electrons + [self.nucleon_pos])

    def update(self, dt, time):
        self.time = time
        if self.balls is not None:
            self.balls.set(pos=self._ball_positions(time))

    def release(self):
        for v in self.clouds.values():
            for c in (v if isinstance(v, list) else [v]):
                c.release()
        for line in self.lines:
            line.release()
        if self.balls is not None:
            self.balls.release()
        self.clouds, self.lines, self.balls, self.electrons = {}, [], None, []

    # ------------------------------------------------------------ dibujo
    def render(self, renderer, cam, didactic, fade):
        if self.mode == "atomic":
            return self._render_atomic(renderer, cam, fade)
        cloud_style = PointStyle(twinkle=0.5, min_px=1.5, max_px=5, softness=0.8)
        cloud_fade = fade * (PLASTIC_CLOUD_FADE if self.balls is not None else 1.0)
        renderer.points(self.clouds["main"], cam, cloud_style, cloud_fade)
        for c in self.clouds.get("neighbors", []):
            renderer.points(c, cam, cloud_style, cloud_fade)
        if "background" in self.clouds:
            renderer.points(self.clouds["background"], cam,
                            PointStyle(min_px=1.5, max_px=300, softness=1.0), fade * 0.6)
        for line in self.lines:
            renderer.lines(line, cam, fade, 0.0)
        # El núcleo: en modo realista, tamaño real (invisible hasta acercarse mucho).
        # En modo didáctico, un marcador mínimo de 5 px para saber dónde está.
        if self.balls is not None:
            self._render_balls(renderer, cam, fade)
            return
        renderer.points(self.clouds["nucleus"], cam,
                        PointStyle(min_px=5.0 if didactic else 0.5, max_px=900, softness=0.6), fade)

    def _render_balls(self, renderer, cam, fade):
        renderer.sprites(self.balls, cam, fade)

    def _render_atomic(self, renderer, cam, fade):
        style = PointStyle(twinkle=0.5, min_px=1.5, max_px=7, softness=0.85)
        main = self.clouds["detail"] if self.view is not None else self.clouds["main"]
        renderer.points(main, cam, style, fade)
        # con un subnivel elegido, los vecinos casi desaparecen: se ve solo lo elegido
        neighbor_fade = fade * (0.15 if self.view is not None else 1.0)
        for c in self.clouds.get("neighbors", []):
            renderer.points(c, cam, style, neighbor_fade)
        if "background" in self.clouds:
            renderer.points(self.clouds["background"], cam,
                            PointStyle(min_px=1.5, max_px=300, softness=1.0), neighbor_fade * 0.6)
        if "axes" in self.clouds:
            renderer.lines(self.clouds["axes"], cam, fade * 0.8, 0.0)
        for line in self.lines:
            renderer.lines(line, cam, neighbor_fade, 0.0)
        renderer.points(self.clouds["nucleus"], cam,
                        PointStyle(min_px=5.0, max_px=900, softness=0.6), fade)

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

    def exit_radius(self) -> float:
        """Distancia al núcleo (pm) a la que se vuelve a la hélice. En el modo atómico es
        mayor: así se puede retroceder (S) y ver el átomo entero, como en AtomicWorld."""
        fraction = ATOMIC_EXIT_FRACTION if self.mode == "atomic" else EXIT_FRACTION
        return fraction * self.vdw_pm

    def ascend_point(self, pos):
        r = np.linalg.norm(pos)
        direction = pos / r if r > 1e-9 else np.array([1.0, 0.0, 0.0])
        return direction * self.exit_radius() * 1.5

    def should_exit(self, pos):
        r = np.linalg.norm(pos)
        if r > self.exit_radius():
            return True
        # si un vecino está más cerca que nuestro núcleo, salimos para entrar en él
        return any(np.linalg.norm(pos - self.rel_pm[j]) < r for j in self.bonded)

    def focus(self, pos, forward, didactic):
        tr, t = self.tr, self.type
        element = tr.element(t.number)
        r = float(np.linalg.norm(pos))
        to_nuc = -pos / max(r, 1e-12)
        if self.mode == "atomic" and self.view is not None and r > 5:
            return self._subshell_focus(r)        # con un subnivel elegido (tecla O), manda él
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
            explanation = {"atomic": "microscope.atom.atomic.cloud_explanation",
                           "plastic": "microscope.atom.plastic.explanation"}.get(
                self.mode, "microscope.atom.cloud_explanation")
            expl = ([tr("microscope.atom.config", config=config_string(t.shells), zeff=zeff)]
                    + tr.lines(explanation))
        else:
            expl = tr.lines("microscope.atom.cloud_explanation_realistic")
        return Focus(tr("microscope.atom.cloud_title", electrons=tr.count(t.electrons, "electron")),
                     atom_context(tr, self.mol, self.i), expl,
                     pair_index=int(self.mol.pair_index[self.i]), distance=r)

    def _subshell_focus(self, r):
        """Modo atómico, tecla O: el subnivel (o el orbital, tecla I) que se está viendo."""
        tr, sh = self.tr, self.shells[self.view]
        zeff = self._zeff(sh)
        if sh.l == 0:
            boxes = "s " + _boxes(sh.electrons)
        else:
            names = orbitals.orbital_names(sh.l)
            occ = orbitals.hund_occupancy(sh.l, sh.electrons)
            boxes = " · ".join((f"[{n} {_boxes(e)}]" if k == self.orbital else f"{n} {_boxes(e)}")
                               for k, (n, e) in enumerate(zip(names, occ)))
        expl = [tr("microscope.atom.atomic.subshell_zeff", zeff=tr.decimal(zeff, 2),
                   r=tr.decimal(orbitals.mean_radius_pm(sh.n, sh.l, zeff), 0)),
                tr(f"microscope.atom.atomic.shape.{'spdf'[sh.l]}", nodes=sh.n - sh.l - 1)]
        if sh.l > 0:
            expl += [tr("microscope.atom.atomic.hund"),
                     tr("microscope.atom.atomic.full" if sh.electrons == 2 * (2 * sh.l + 1)
                        else "microscope.atom.atomic.partial"),
                     tr("microscope.atom.atomic.orbital_colors")]
        if self.orbital is not None:
            e = orbitals.hund_occupancy(sh.l, sh.electrons)[self.orbital]
            title = tr("microscope.atom.atomic.orbital_title", name=self.orbital_label(),
                       electrons=tr.count(e, "electron"))
        else:
            title = tr("microscope.atom.atomic.subshell_title", name=self.subshell_label(self.view),
                       electrons=tr.count(sh.electrons, "electron"))
        return Focus(title, boxes, expl, pair_index=int(self.mol.pair_index[self.i]), distance=r)

    def time_info(self):
        return self.tr("microscope.atom.time")


def _boxes(electrons: int) -> str:
    """Casilla de un orbital: '↑↓', '↑' o '—' (vacío)."""
    return {2: "↑↓", 1: "↑"}.get(electrons, "—")
