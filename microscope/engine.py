"""
Motor: pila de niveles, transiciones con fundido y dibujo.

    pila = [HelixLevel, AtomLevel, NucleusLevel, NucleonLevel]   (el último es donde está la nave)
           [ProteinLevel, AtomLevel, ...]                         (una proteína)

La nave siempre vive en las coordenadas del nivel superior de la pila.
Al entrar en un hijo se apila; al salir se desapila. Durante ~0,6 s se
dibujan ambos niveles con un fundido cruzado.

Comparar dos proteínas (la natural y la tuya, ya superpuestas: comparten coordenadas).
Cada "vista" tiene su propia pila; la nave es la misma para todas:

    split     pantalla partida: natural a la izquierda, la tuya a la derecha, con la misma
              cámara. Al entrar en un átomo, se ve a pantalla completa el lado que manejas
              (tecla V para cambiarlo) y al volver a la proteína se parte otra vez.
    overlay   las dos en la misma escena, cada una de un color (tecla C para alternar).
"""
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np

from microscope.autopilot import Autopilot
from microscope.levels.base import Level
from microscope.levels.helix import HelixLevel
from microscope.mathutil import perspective, view_matrix
from microscope.render.points import Camera, Renderer
from microscope.ship import Ship

FOV_Y = 70.0
TRANSITION_TIME = 0.6
CACHE_SIZE = 4


@dataclass
class Transition:
    old: Level
    to_old: Callable[[np.ndarray], np.ndarray]   # posición actual -> coords del nivel viejo
    t: float = 0.0


@dataclass
class View:
    """Una pila de niveles con su caché y su transición."""
    stack: list                                   # [(nivel, Child por el que se entró)]
    label: str = ""
    cache: OrderedDict = field(default_factory=OrderedDict)
    transition: Optional[Transition] = None

    @property
    def root(self) -> Level:
        return self.stack[0][0]


@dataclass
class Comparison:
    """Lo necesario para comparar dos proteínas en el microscopio."""
    second: object                                # ProteinMolecule (la tuya)
    labels: tuple[str, str]                       # ("Hemoglobin subunit beta · P68871", ...)
    highlight: tuple[dict, dict] = ({}, {})       # aminoácido -> "Glu7Val", en cada una


def make_root(molecule, tr, label="", highlight=None) -> Level:
    if getattr(molecule, "kind", "dna") == "protein":
        from microscope.levels.protein import ProteinLevel
        return ProteinLevel(molecule, tr, label, highlight)
    return HelixLevel(molecule, tr)


class Engine:
    def __init__(self, ctx, molecule, tr, didactic=True, comparison: Comparison | None = None,
                 label: str = "", highlight: dict | None = None):
        self.renderer = Renderer(ctx)
        self.molecule = molecule
        self.tr = tr
        self.didactic = didactic
        self.time = 0.0
        self.comparison = comparison

        if comparison is None:
            root = make_root(molecule, tr, label, highlight)
            root.build(self.renderer, didactic)
            self.views = {"single": View([(root, None)], label)}
            self.layout, self.active = "single", "single"
        else:
            self.views = {}
            sides = (("left", molecule, comparison.labels[0], comparison.highlight[0]),
                     ("right", comparison.second, comparison.labels[1], comparison.highlight[1]))
            for side, mol, text, marks in sides:
                root = make_root(mol, tr, text, marks)
                root.build(self.renderer, didactic)
                self.views[side] = View([(root, None)], text)
            self.layout, self.active = "split", "left"

        # la nave empieza lejos, como si llegara desde el espacio (una proteína, compacta, se
        # ve mejor un poco más cerca: así llena la vista, también en media pantalla)
        extent = np.ptp(molecule.positions, axis=0).max()
        if getattr(molecule, "kind", "dna") == "protein":
            start = np.array([extent * 1.3 + 2.0, -extent * 0.3, extent * 0.2])
        else:
            start = np.array([extent * 2.2 + 6.0, -extent * 0.4, extent * 0.25])
        self.ship = Ship(start, look_at=np.zeros(3))
        self.home = start.copy()                     # vista panorámica (Retroceso en la raíz)
        self.autopilot = Autopilot(tr)

    # ------------------------------------------------------------ vistas
    @property
    def view(self) -> View:
        """La vista que se maneja: la del lado activo, o la superpuesta."""
        return self.views["overlay" if self.layout == "overlay" else self.active]

    @property
    def stack(self) -> list:
        return self.view.stack

    @property
    def cache(self) -> OrderedDict:
        return self.view.cache

    @property
    def transition(self) -> Optional[Transition]:
        return self.view.transition

    @transition.setter
    def transition(self, value: Optional[Transition]):
        self.view.transition = value

    @property
    def at_root(self) -> bool:
        return len(self.stack) == 1

    @property
    def split_visible(self) -> bool:
        """Pantalla partida a la vista: modo split y en el nivel de la proteína."""
        return self.layout == "split" and self.at_root

    def toggle_layout(self) -> bool:
        """Partida <-> superpuestas (solo desde el nivel de la proteína). False si no se puede."""
        if self.comparison is None or not self.at_root:
            return False
        if self.layout == "split":
            if "overlay" not in self.views:
                from microscope.levels.protein import OverlayLevel, ProteinLevel
                c = self.comparison
                levels = [ProteinLevel(self.molecule, self.tr, c.labels[0], c.highlight[0]),
                          ProteinLevel(c.second, self.tr, c.labels[1], c.highlight[1])]
                overlay = OverlayLevel(levels, self.tr)
                overlay.build(self.renderer, self.didactic)
                self.views["overlay"] = View([(overlay, None)])
            self.layout = "overlay"
        else:
            self.layout = "split"
        self.transition = None
        return True

    def switch_side(self) -> bool:
        """Cambia el lado que manejas (pantalla partida, en el nivel de la proteína)."""
        if not self.split_visible:
            return False
        self.active = "right" if self.active == "left" else "left"
        return True

    # ------------------------------------------------------------ propiedades
    @property
    def level(self) -> Level:
        return self.stack[-1][0]

    def surface_distance(self) -> float:
        return self.level.surface_distance(self.ship.pos)

    # ------------------------------------------------------------ lógica
    def update(self, dt, move, boost):
        self.time += dt
        self.level.update(dt, self.time)
        if self.transition:
            self.transition.old.update(dt, self.time)
        factor = 1.0
        if self.autopilot.active:
            move, factor = self.autopilot.step(self, dt)
            boost = False
            if self.autopilot.mode in ("face", None):
                self.ship.vel *= np.exp(-dt * 8.0)      # frenar al girarse
        self.ship.update(dt, move, boost, self.surface_distance(), factor)
        self._check_level_change()
        if self.transition:
            self.transition.t += dt / TRANSITION_TIME
            if self.transition.t >= 1:
                self.transition = None

    def _check_level_change(self):
        level = self.level
        if len(self.stack) > 1 and level.should_exit(self.ship.pos):
            child = self.stack[-1][1]
            self.stack.pop()
            self.ship.unscale(child.center, child.scale)
            self.transition = Transition(level, lambda p, c=child: (p - c.center) * c.scale)
            self.autopilot.on_level_change(self, entered=False, exited_center=child.center)
            return
        only = self.autopilot.only_key
        if only == "__none__":           # el piloto va a un destino sin hijo (un quark)
            return
        child = level.child_to_enter(self.ship.pos, only)
        if child is not None:
            new = self._get_level(child)
            self.ship.rescale(child.center, child.scale)
            self.stack.append((new, child))
            self.transition = Transition(level, lambda p, c=child: p / c.scale + c.center)
            self.autopilot.on_level_change(self, entered=True)

    def _get_level(self, child) -> Level:
        if child.key in self.cache:
            self.cache.move_to_end(child.key)
            return self.cache[child.key]
        lvl = child.create()
        lvl.build(self.renderer, self.didactic)
        self.cache[child.key] = lvl
        while len(self.cache) > CACHE_SIZE:
            key, old = next(iter(self.cache.items()))
            if any(old is l for l, _ in self.stack):
                self.cache.move_to_end(key)
                break
            self.cache.popitem(last=False)
            old.release()
        return lvl

    def toggle_mode(self):
        self.didactic = not self.didactic
        for view in self.views.values():
            levels = [l for l, _ in view.stack] + [l for l in view.cache.values()
                                                     if not any(l is s for s, _ in view.stack)]
            for lvl in levels:
                lvl.set_mode(self.renderer, self.didactic)
            view.transition = None

    # ------------------------------------------------------------ dibujo
    def _camera(self, level, pos, width, height) -> Camera:
        d = max(level.surface_distance(pos), 1e-12)
        near = d * 0.02
        far = near * 1e7
        proj = perspective(FOV_Y, width / max(height, 1), near, far)
        return Camera(view_matrix(pos, self.ship.orient), proj, height, self.time)

    def render(self, width, height):
        r = self.renderer
        r.ctx.viewport = (0, 0, width, height)
        r.begin()
        if self.split_visible:
            half = width // 2
            for x, side in ((0, "left"), (half, "right")):
                self._render_view(self.views[side], x, half, height)
            r.ctx.viewport = (0, 0, width, height)
        else:
            self._render_view(self.view, 0, width, height)

    def _render_view(self, view: View, x: int, width: int, height: int):
        r = self.renderer
        r.ctx.viewport = (x, 0, width, height)
        level = view.stack[-1][0]
        fade_new = 1.0
        if view.transition:
            tr = view.transition
            s = min(tr.t, 1.0)
            fade_new = s * s * (3 - 2 * s)     # suavizado
            old_pos = tr.to_old(self.ship.pos)
            tr.old.render(r, self._camera(tr.old, old_pos, width, height),
                          self.didactic, 1.0 - fade_new)
        level.render(r, self._camera(level, self.ship.pos, width, height), self.didactic,
                     fade_new)

    # ------------------------------------------------------------ datos para el HUD
    def hud_state(self, width, height) -> dict:
        lvl = self.level
        d = self.surface_distance()
        aspect = width / max(height, 1)
        visible_w = 2 * d * np.tan(np.radians(FOV_Y) / 2) * aspect * lvl.unit_m
        speed = np.linalg.norm(self.ship.vel) * lvl.unit_m
        return {
            "level": lvl,
            "depth_m": d * lvl.unit_m,
            "visible_width_m": visible_w,
            "speed_m_s": speed,
            "focus": lvl.focus(self.ship.pos, self.ship.forward, self.didactic),
            "time": lvl.time_info(),
            "didactic": self.didactic,
            "molecule": self.view.root.mol if hasattr(self.view.root, "mol") else self.molecule,
            "panels": self._panels(),
            "comparing": self.comparison is not None,
            "speed_mult": self.ship.speed_mult,
            "autopilot": self.autopilot.status(),
            "target_screen": self._project(self.autopilot.target_point(), width, height),
        }

    def _panels(self) -> list[dict]:
        """Etiquetas de cada mitad para el HUD ([] si no hay pantalla partida a la vista).
        x0, x1: fracción de la anchura que ocupa cada mitad."""
        if not self.split_visible:
            return []
        return [{"label": self.views[side].label, "x0": x0, "x1": x0 + 0.5,
                 "active": side == self.active}
                for side, x0 in (("left", 0.0), ("right", 0.5))]

    def _project(self, point, width, height):
        """Posición en pantalla (píxeles, origen arriba-izquierda) o None si está detrás."""
        if point is None or not self.autopilot.active:
            return None
        x0, w = 0, width
        if self.split_visible:               # el destino se marca en la mitad que manejas
            w = width // 2
            x0 = 0 if self.active == "left" else w
        cam = self._camera(self.level, self.ship.pos, w, height)
        clip = cam.proj @ cam.view @ np.append(point, 1.0)
        if clip[3] <= 0:
            return None
        ndc = clip[:3] / clip[3]
        return (x0 + (ndc[0] * 0.5 + 0.5) * w, (0.5 - ndc[1] * 0.5) * height)
