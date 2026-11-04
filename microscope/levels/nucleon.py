"""
Nivel 4: dentro de un protón o un neutrón (unidad: 0,1 fm = 10⁻¹⁶ m).

Física (cromodinámica cuántica, QCD):
  - 3 quarks de valencia: protón = u u d, neutrón = u d d. Sus masas salen de
    atomicworld/data/particles.json (PDG 2024).
  - Mar de partones: pares quark–antiquark y gluones que aparecen y desaparecen sin parar.
  - Carga de color: cada quark es rojo, verde o azul (etiquetas, no colores reales). El
    nucleón siempre es "blanco". Los quarks intercambian color emitiendo y absorbiendo gluones.
  - Confinamiento: los quarks están unidos por un tubo de flujo de gluones en forma de Y
    (resultado de la QCD en retículo). Nunca se ha observado un quark aislado.
  - Los quarks son puntuales: tamaño < 10⁻¹⁸ m. Aquí termina el viaje.

Modos:
  - Realista: solo la "niebla" de partones con el perfil de carga medido (radio 0,84 fm).
  - Didáctico: los 3 quarks con su color, el tubo de gluones en Y y los gluones viajando
    entre quarks e intercambiando el color.

El movimiento de los quarks y el tubo en Y son una ilustración: la QCD no tiene trayectorias.
"""
from functools import lru_cache

import numpy as np

from atomicworld.repository import ParticleCatalog
from microscope.levels.base import Focus, Level, Target, aimed
from microscope.render.points import PointCloud, PointStyle

UNIT_PER_FM = 10.0
NUCLEON_A = 0.84 / np.sqrt(12) * UNIT_PER_FM     # perfil exponencial (rms 0,84 fm)
QUARK_ORBIT = 3.2                                 # amplitud del movimiento (unidades)
EXIT_RADIUS = 13.0                                # se sale a 1,3 fm del centro
QUARK_STOP = 2.5                                  # el piloto automático se para a 0,25 fm del quark

COLORS = {"r": (1.00, 0.22, 0.22), "g": (0.25, 1.00, 0.35), "b": (0.30, 0.45, 1.00)}
ANTI = {"r": (0.20, 1.00, 1.00), "g": (1.00, 0.30, 1.00), "b": (1.00, 1.00, 0.25)}
REALISTIC_COLOR = (0.95, 0.82, 0.72)
GLUON_COLOR = (0.85, 0.85, 0.95)
CONTENT = {"p": "uud", "n": "udd"}
FLAVOR_NAMES = {"u": "Up", "d": "Down"}

SEA_POINTS = 140_000
TUBE_POINTS = 1500       # por brazo de la Y
GLUON_POINTS = 400       # gluón viajero
EXCHANGE_EVERY = 1.3     # s entre intercambios de color
EXCHANGE_TIME = 0.55     # s que tarda el gluón en viajar


@lru_cache(maxsize=1)
def quark_masses_mev() -> dict:
    catalog = ParticleCatalog()
    return {flavor: catalog.find(name).mass_mev for flavor, name in FLAVOR_NAMES.items()}


def _directions(n, rng):
    v = rng.normal(size=(n, 3))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


class NucleonLevel(Level):
    depth = 3
    unit_m = 1e-16

    def __init__(self, kind: str, seed: int, tr, neighbors=None, neighbor_kinds=None):
        super().__init__(tr)
        self.kind = kind
        self.flavors = list(CONTENT[kind])
        self.rng = np.random.default_rng(seed)
        self.colors = list("rgb")                     # siempre uno de cada: "blanco"
        self.rng.shuffle(self.colors)
        # trayectorias de Lissajous distintas para cada quark
        self.freq = self.rng.uniform(0.6, 1.3, (3, 3))
        self.phase = self.rng.uniform(0, 2 * np.pi, (3, 3))
        self.neighbors = np.zeros((0, 3)) if neighbors is None else neighbors
        self.neighbor_kinds = neighbor_kinds or []
        self.quark_pos = self._quarks(0.0)
        self.exchange = None                          # (i, j, t0)
        self.next_exchange = EXCHANGE_EVERY
        self.didactic = True
        self.time = 0.0
        self.static = []
        self.dynamic = None
        # desplazamientos fijos del tubo y del gluón (se recolocan cada frame)
        self.tube_t = self.rng.random(3 * TUBE_POINTS)
        self.tube_off = self.rng.normal(0, 0.3, (3 * TUBE_POINTS, 3))
        self.gluon_off = self.rng.normal(0, 0.35, (GLUON_POINTS, 3))

    @property
    def name(self):
        return self.tr("microscope.levels." + ("proton" if self.kind == "p" else "neutron"))

    # ------------------------------------------------------------ física de la animación
    def _quarks(self, t):
        q = QUARK_ORBIT * np.sin(self.freq * t + self.phase)
        return q - q.mean(axis=0) * 0.8           # el centro de masas casi quieto

    def update(self, dt, time):
        self.time = time
        self.quark_pos = self._quarks(time)
        if self.exchange is None and time >= self.next_exchange:
            i, j = self.rng.choice(3, 2, replace=False)
            self.exchange = (int(i), int(j), time)
        if self.exchange and time - self.exchange[2] >= EXCHANGE_TIME:
            i, j, _ = self.exchange
            # el gluón lleva color de i y anticolor de j: al absorberlo, intercambian color
            self.colors[i], self.colors[j] = self.colors[j], self.colors[i]
            self.exchange = None
            self.next_exchange = time + EXCHANGE_EVERY * self.rng.uniform(0.6, 1.4)
        if self.dynamic is not None and self.didactic:
            self.dynamic.write(*self._dynamic_points())

    # ------------------------------------------------------------ construcción
    def build(self, renderer, didactic):
        self.didactic = didactic
        rng = np.random.default_rng(11)
        n = SEA_POINTS
        r = rng.gamma(3.0, NUCLEON_A, n)
        pos = _directions(n, rng) * r[:, None]
        if didactic:
            col = np.tile(GLUON_COLOR, (n, 1))
            sea = rng.random(n) < 0.35                     # 35 %: quarks/antiquarks del mar
            pal = np.array(list(COLORS.values()) + list(ANTI.values()))
            col[sea] = pal[rng.integers(0, 6, sea.sum())]
            alpha = 0.045
        else:
            col = np.tile(REALISTIC_COLOR, (n, 1))
            alpha = 0.14
        col = np.hstack([col, np.full((n, 1), alpha)])
        self.static.append(("sea", PointCloud(renderer, pos, col, np.full(n, 0.25),
                                              rng.random(n).astype("f4") * 1000)))

        # nucleones vecinos de fondo (a ~1,8 fm)
        if len(self.neighbors):
            pts, cols = [], []
            for c, k in zip(self.neighbors, self.neighbor_kinds):
                m = 2500
                pts.append(c + _directions(m, rng) * rng.gamma(3.0, NUCLEON_A, m)[:, None])
                rgb = ((1.0, 0.36, 0.30) if k == "p" else (0.55, 0.68, 0.90)) if didactic \
                    else REALISTIC_COLOR
                cols.append(np.tile([*rgb, 0.08], (m, 1)))
            pts = np.vstack(pts)
            self.static.append(("neighbors", PointCloud(
                renderer, pts, np.vstack(cols), np.full(len(pts), 0.3),
                rng.random(len(pts)).astype("f4") * 1000)))

        if didactic:
            self.dynamic = PointCloud(renderer, *self._dynamic_points(), dynamic=True)

    def _dynamic_points(self):
        q = self.quark_pos
        junction = q.mean(axis=0)
        pos, col, size = [], [], []
        # tubo de flujo en Y: cada brazo va del quark al punto de unión
        for a in range(3):
            sl = slice(a * TUBE_POINTS, (a + 1) * TUBE_POINTS)
            t = self.tube_t[sl][:, None]
            pos.append(q[a] * (1 - t) + junction * t + self.tube_off[sl] * (0.4 + 0.6 * t))
            c = np.array(COLORS[self.colors[a]]) * 0.5 + 0.5
            col.append(np.tile([*c, 0.22], (TUBE_POINTS, 1)))
            size.append(np.full(TUBE_POINTS, 0.2))
        # gluón viajando de i a j (color de i + anticolor de j)
        if self.exchange:
            i, j, t0 = self.exchange
            s = np.clip((self.time - t0) / EXCHANGE_TIME, 0, 1)
            center = q[i] * (1 - s) + q[j] * s
            half = GLUON_POINTS // 2
            c1 = [*COLORS[self.colors[i]], 0.7]
            c2 = [*ANTI[self.colors[j]], 0.7]
            pos.append(center + self.gluon_off)
            col.append(np.vstack([np.tile(c1, (half, 1)), np.tile(c2, (GLUON_POINTS - half, 1))]))
        else:
            pos.append(np.tile(q.mean(axis=0), (GLUON_POINTS, 1)))
            col.append(np.zeros((GLUON_POINTS, 4)))
        size.append(np.full(GLUON_POINTS, 0.15))
        # los 3 quarks de valencia: puntuales, se dibujan con un brillo mínimo
        pos.append(q)
        col.append(np.array([[*COLORS[c], 1.0] for c in self.colors]))
        size.append(np.full(3, 0.02))
        pos, col, size = np.vstack(pos), np.vstack(col), np.concatenate(size)
        return pos, col, size, np.arange(len(pos), dtype="f4")

    def release(self):
        for _, c in self.static:
            c.release()
        if self.dynamic is not None:
            self.dynamic.release()
        self.static, self.dynamic = [], None

    # ------------------------------------------------------------ dibujo
    def render(self, renderer, cam, didactic, fade):
        for name, cloud in self.static:
            tw = 0.9 if name == "sea" else 0.3
            renderer.points(cloud, cam, PointStyle(twinkle=tw, min_px=1.5, max_px=5,
                                                   softness=0.8), fade)
        if didactic and self.dynamic is not None:
            n = self.dynamic.count - 3
            renderer.points(self.dynamic, cam, PointStyle(min_px=1.5, max_px=5, softness=0.8),
                            fade, first=0, count=n)
            # los quarks no tienen tamaño: brillo mínimo de 12 px para poder verlos
            renderer.points(self.dynamic, cam, PointStyle(min_px=16.0, max_px=20,
                                                          softness=1.0, dim_small=False), fade, first=n, count=3)

    # ------------------------------------------------------------ navegación
    def surface_distance(self, pos):
        d = np.linalg.norm(self.quark_pos - pos, axis=1).min()
        return float(max(min(d, np.linalg.norm(pos) + 1.0), 0.05))

    def should_exit(self, pos):
        return np.linalg.norm(pos) > EXIT_RADIUS

    def dive_target(self, pos, forward, didactic):
        """No hay nivel más abajo: el piloto se acerca al quark (que se mueve) y se para."""
        k = aimed(self.quark_pos, pos, forward, 6)
        if k is None:
            k = int(np.argmin(np.linalg.norm(self.quark_pos - pos, axis=1)))
        flavor = "up" if self.flavors[k] == "u" else "down"
        return Target(lambda: self.tr(f"microscope.nucleon.{flavor}.title"),
                      lambda: self.quark_pos[k], key=None, arrive_radius=QUARK_STOP)

    def ascend_point(self, pos):
        r = np.linalg.norm(pos)
        direction = pos / r if r > 1e-9 else np.array([1.0, 0.0, 0.0])
        return direction * EXIT_RADIUS * 1.5

    def focus(self, pos, forward, didactic):
        tr = self.tr
        name = tr("microscope.nucleon." + ("proton" if self.kind == "p" else "neutron"))
        content = CONTENT[self.kind]
        color_name = lambda c: tr(f"microscope.nucleon.colors.{c}")
        if didactic:
            v = self.quark_pos - pos
            d = np.linalg.norm(v, axis=1)
            cos = (v @ forward) / np.maximum(d, 1e-9)
            score = np.where((cos > np.cos(np.radians(5))) | (d < QUARK_STOP * 1.4), d, np.inf)
            k = int(np.argmin(score))
            if np.isfinite(score[k]):
                flavor = "up" if self.flavors[k] == "u" else "down"
                mass = quark_masses_mev()[self.flavors[k]]
                return Focus(tr(f"microscope.nucleon.{flavor}.title"),
                             tr(f"microscope.nucleon.{flavor}.detail", mass=tr.decimal(mass, 2),
                                  color=color_name(self.colors[k])),
                             tr.lines("microscope.nucleon.quark_explanation"), distance=float(d[k]))
            expl = tr.lines("microscope.nucleon.explanation")
            detail = tr("microscope.nucleon.colors_now", colors=", ".join(color_name(c) for c in self.colors))
        else:
            expl = tr.lines("microscope.nucleon.explanation_realistic")
            detail = tr("microscope.nucleon.realistic_detail")
        return Focus(tr("microscope.nucleon.title", name=name, content=" ".join(content)),
                     detail, expl, distance=float(np.linalg.norm(pos)))

    def time_info(self):
        return self.tr("microscope.nucleon.time")
