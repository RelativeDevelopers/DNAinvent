"""
Nivel 3: el núcleo atómico (unidad: femtómetro, 1 fm = 10⁻¹⁵ m).

Modos:
  - Didáctico: protones (rojo) y neutrones (gris azulado) por separado, con líneas para la
    fuerza nuclear fuerte entre vecinos (alcance ~2 fm).
  - Realista: no se puede decir "dónde" está cada nucleón (son ondas cuánticas). Se dibuja
    la densidad del núcleo entero con el perfil de Woods–Saxon. Para el hidrógeno-1 (un
    único protón) se dibuja ese protón.

  - Plástico: cada protón (+, rojo) y cada neutrón (0, azul) es una bola de plástico
    (render/sprites.py), con un pequeño temblor. Al acercarte a una se vuelve transparente
    para poder entrar.

La física (posiciones, perfiles) es la de AtomicWorld: atomicworld/scene/nucleus.py. Con la
misma semilla, el núcleo del carbono-12 es idéntico en las dos apps.
"""
import numpy as np

from atomicworld import atoms
from atomicworld.constants import NEUTRON_MASS_MEV, PROTON_MASS_MEV
from atomicworld.scene import nucleus as nuc
from atomicworld.scene.orbitals import probability_inside_nucleus
from microscope.levels import plastic
from microscope.levels.base import Child, Focus, Level, Target, aimed
from microscope.render.points import LineSet, PointCloud, PointStyle

PROTON_COLOR = (1.00, 0.36, 0.30)
NEUTRON_COLOR = (0.55, 0.68, 0.90)
REALISTIC_COLOR = (0.95, 0.82, 0.72)
REALISTIC_POINTS = 160_000
VIB_AMP_FM = 0.22          # movimiento de Fermi de los nucleones (ver time_info)
ENTER_NUCLEON_FM = 0.45    # se entra en un nucleón a 0,45 fm de su centro
EXIT_FACTOR = 4.0          # se sale al átomo a 4 × R nuclear
UNIT_PER_FM = 10.0         # el nivel del nucleón trabaja en décimas de fm
NEAR_FADE = 2.0            # modo plástico: la bola se vuelve transparente a menos de 2 radios


class NucleusLevel(Level):
    depth = 2
    unit_m = 1e-15

    def __init__(self, atom_type, tr, context=lambda: ""):
        super().__init__(tr)
        self.type = atom_type
        self.el, self.iso = atom_type.element, atom_type.isotope
        self.context = context           # 'Adenina (A) · par 3 · hebra 1' (función: se traduce)
        self.A = self.iso.a
        self.R = atoms.nuclear_radius_fm(self.A)
        rng = np.random.default_rng(self.el.number * 1000 + self.A)     # igual que AtomicWorld
        self.centers = nuc.pack_nucleons(self.A, self.R, rng)
        self.kinds = nuc.assign_kinds(self.el.protons, self.iso.neutrons, rng)
        self.clouds = []
        self.lines = []
        self.balls = None                # modo plástico (SpriteSet)
        self.wobble = plastic.Wobble(self.A, plastic.FERMI_WOBBLE_FM, np.random.default_rng(self.A))

    @property
    def name(self):
        tr = self.tr
        element = tr.element(self.el.number)
        return tr("microscope.levels.nucleus", element=element, of_element=tr.of(element), a=self.A)

    # ------------------------------------------------------------ construcción
    def build(self, renderer, didactic):
        if self.mode == "plastic":
            from microscope.render.sprites import SpriteSet
            self.balls = SpriteSet(renderer, self.centers, np.full(self.A, plastic.NUCLEON_BALL_FM),
                                   plastic.nucleon_kinds(self.kinds), min_px=2.0)
            return
        rng = np.random.default_rng(7)
        if didactic or self.A == 1:
            pos, owner = nuc.nucleon_clouds(self.centers, nuc.points_per_nucleon(self.A), rng)
            if didactic:
                rgb = np.where((self.kinds[owner] == "p")[:, None], PROTON_COLOR, NEUTRON_COLOR)
            else:
                rgb = np.tile(REALISTIC_COLOR, (len(pos), 1))
            col = np.hstack([rgb, np.full((len(pos), 1), 0.22)])
            # misma semilla en todo el nucleón: vibra como un bloque
            self.clouds.append(("nucleons", PointCloud(renderer, pos, col, np.full(len(pos), 0.03),
                                                       owner.astype("f4"))))
        else:
            n = REALISTIC_POINTS
            r = nuc.sample_woods_saxon(n, self.R, nuc.WS_DIFFUSENESS_FM, rng)
            pts = nuc.directions(n, rng) * r[:, None]
            col = np.tile([*REALISTIC_COLOR, 0.18], (n, 1))
            self.clouds.append(("matter", PointCloud(renderer, pts, col, np.full(n, 0.03),
                                                     rng.random(n).astype("f4") * 1000)))
        if didactic and self.A > 1:
            a, b = nuc.strong_bonds(self.centers)
            self.lines.append(LineSet(renderer, a, b, np.tile([0.45, 1.0, 0.55, 0.75], (len(a), 1))))

    def release(self):
        for _, c in self.clouds:
            c.release()
        for line in self.lines:
            line.release()
        if self.balls is not None:
            self.balls.release()
        self.clouds, self.lines, self.balls = [], [], None

    def update(self, dt, time):
        if self.balls is not None:
            self.balls.set(pos=self.centers + self.wobble.offsets(time))

    # ------------------------------------------------------------ dibujo
    def render(self, renderer, cam, didactic, fade):
        if self.balls is not None:
            renderer.sprites(self.balls, cam, fade, near_fade=NEAR_FADE)
        for kind, cloud in self.clouds:
            if kind == "nucleons":
                style = PointStyle(vib_amp=VIB_AMP_FM, min_px=1.5, max_px=5, softness=0.8)
            else:
                style = PointStyle(twinkle=0.4, min_px=1.5, max_px=5, softness=0.8)
            renderer.points(cloud, cam, style, fade)
        for line in self.lines:
            renderer.lines(line, cam, fade, 0.12)

    # ------------------------------------------------------------ navegación
    def surface_distance(self, pos):
        d = np.linalg.norm(self.centers - pos, axis=1).min()
        return float(max(d * 0.6, 0.05))

    def child_to_enter(self, pos, only=None):
        if only is not None and only[1] != id(self):
            return None
        d = np.linalg.norm(self.centers - pos, axis=1)
        k = int(np.argmin(d)) if only is None else only[2]
        if d[k] >= ENTER_NUCLEON_FM:
            return None
        from microscope.levels.nucleon import NucleonLevel
        others = [j for j in range(self.A) if j != k]
        rel = (self.centers[others] - self.centers[k]) * UNIT_PER_FM
        kinds = [str(self.kinds[j]) for j in others]
        kind = str(self.kinds[k])
        return Child(key=("nucleon", id(self), k), center=self.centers[k].copy(), scale=UNIT_PER_FM,
                     create=lambda: NucleonLevel(kind, k, self.tr, neighbors=rel, neighbor_kinds=kinds))

    def dive_target(self, pos, forward, didactic):
        k = aimed(self.centers, pos, forward, 6)
        if k is None:      # ninguno delante: el más cercano a la nave
            k = int(np.argmin(np.linalg.norm(self.centers - pos, axis=1)))
        p = self.centers[k].copy()
        kind = "proton" if self.kinds[k] == "p" else "neutron"
        return Target(lambda: self.tr(f"microscope.nucleus.{kind}.title"), lambda: p,
                      key=("nucleon", id(self), k))

    def ascend_point(self, pos):
        r = np.linalg.norm(pos)
        direction = pos / r if r > 1e-9 else np.array([1.0, 0.0, 0.0])
        return direction * EXIT_FACTOR * self.R * 1.5

    def should_exit(self, pos):
        return np.linalg.norm(pos) > EXIT_FACTOR * self.R

    def focus(self, pos, forward, didactic):
        tr = self.tr
        k = aimed(self.centers, pos, forward, 6, close=nuc.NUCLEON_RMS_FM)
        realistic_blob = not didactic and self.A > 1
        if k is not None and not realistic_blob:
            kind = "proton" if self.kinds[k] == "p" else "neutron"
            mass = PROTON_MASS_MEV if kind == "proton" else NEUTRON_MASS_MEV
            expl = []
            if didactic:
                expl = tr.lines(f"microscope.nucleus.{kind}.explanation")
                if self.mode == "plastic":
                    expl += tr.lines("microscope.nucleus.plastic")
                elif self.A > 1:
                    expl.append(tr("microscope.nucleus.strong_lines"))
            detail = tr(f"microscope.nucleus.{kind}.detail", mass=tr.decimal(mass, 2))
            return Focus(tr(f"microscope.nucleus.{kind}.title"),
                         detail + " · " + tr("microscope.nucleus.go_in"), expl,
                         distance=float(np.linalg.norm(self.centers[k] - pos)))

        p_in = probability_inside_nucleus(self.el.protons, self.type.shells, self.R / 1000)
        protons = tr.count(self.el.protons, "proton")
        expl = []
        if didactic:
            key = "explanation_hydrogen" if self.el.protons == 1 else "explanation"
            expl = tr.lines(f"microscope.nucleus.{key}", protons=protons, p_in=tr.sci(p_in))
            if self.mode == "plastic":
                expl += tr.lines("microscope.nucleus.plastic")
        elif self.A > 1:
            expl = tr.lines("microscope.nucleus.explanation_realistic")
        title = tr("microscope.nucleus.title", protons=protons,
                   neutrons=tr.count(self.iso.neutrons, "neutron"), r=tr.decimal(self.R, 1))
        detail = tr("microscope.nucleus.detail", rho=tr.sci(nuc.density_kg_m3(self.A, self.R)))
        context = self.context()
        return Focus(title, f"{context} · {detail}" if context else detail, expl,
                     distance=float(np.linalg.norm(pos)))

    def time_info(self):
        return self.tr("microscope.nucleus.time")
