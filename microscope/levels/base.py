"""
Interfaz común de todos los niveles de profundidad.

Cada nivel trabaja en su propia unidad (nm, pm, fm, 0,1 fm) para que los números siempre
estén en un rango cómodo. Al entrar en un "hijo", el motor convierte la posición de la
nave a las coordenadas locales del hijo:

    pos_hijo = (pos_padre - child.center) * child.scale

Los textos se piden al traductor (self.tr) en el momento de mostrarlos: así la tecla L
cambia el idioma al instante, sin reconstruir nada.
"""
from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np


@dataclass
class Child:
    key: object                      # identificador estable (para cachear)
    center: np.ndarray               # en unidades del padre
    scale: float                     # unidades del hijo por unidad del padre
    create: Callable[[], "Level"]    # construye el nivel interior


@dataclass
class Target:
    """Destino del piloto automático dentro de un nivel."""
    label: Callable[[], str]                     # texto para el HUD (se traduce al pedirlo)
    point: Callable[[], np.ndarray]              # posición actual (puede moverse: quarks)
    key: object = None                           # Child.key en el que entrar (None = no hay hijo)
    arrive_radius: float = 0.0                   # sin hijo: distancia a la que se considera llegada


@dataclass
class Focus:
    """Lo que la nave tiene delante, para el HUD."""
    title: str
    detail: str = ""
    explanation: list = field(default_factory=list)   # solo modo didáctico
    pair_index: Optional[int] = None                   # para resaltar en la secuencia
    distance: Optional[float] = None                   # unidades locales


class Level:
    depth = 0             # 0 = hélice, 1 = átomo, 2 = núcleo, 3 = nucleón
    unit_m = 1.0          # metros por unidad local
    # modo de visualización (engine.MODES). Lo pone el motor antes de build(); el argumento
    # didactic de build/render/focus sigue diciendo si los textos son didácticos.
    mode = "didactic"

    def __init__(self, tr):
        self.tr = tr

    @property
    def name(self) -> str:
        raise NotImplementedError

    def build(self, renderer, didactic: bool):
        raise NotImplementedError

    def set_mode(self, renderer, didactic: bool):
        self.release()
        self.build(renderer, didactic)

    def update(self, dt: float, time: float):
        """Animación propia del nivel (opcional)."""

    def render(self, renderer, cam, didactic: bool, fade: float):
        raise NotImplementedError

    def surface_distance(self, pos: np.ndarray) -> float:
        """Distancia (unidades locales) al objeto más cercano: fija la velocidad de la nave."""
        raise NotImplementedError

    def child_to_enter(self, pos: np.ndarray, only=None) -> Optional[Child]:
        """only: si se indica, solo se puede entrar en el hijo con esa clave
        (el piloto automático no debe 'caer' en otro átomo por el camino)."""
        return None

    def dive_target(self, pos, forward, didactic) -> Optional[Target]:
        """Piloto automático (Enter): a qué hijo ir. Prioriza lo que hay en la mira."""
        return None

    def ascend_point(self, pos) -> Optional[np.ndarray]:
        """Piloto automático (Retroceso): un punto fuera del nivel por donde salir.
        None = nivel raíz (el motor vuelve a la vista inicial)."""
        return None

    def should_exit(self, pos: np.ndarray) -> bool:
        return False

    def focus(self, pos: np.ndarray, forward: np.ndarray, didactic: bool) -> Optional[Focus]:
        return None

    def time_info(self) -> str:
        return ""

    def release(self):
        pass


def aimed(points: np.ndarray, pos: np.ndarray, forward: np.ndarray, cone_deg: float,
          close: float = 0.0) -> Optional[int]:
    """Índice del punto más cercano dentro del cono de la mira (o a menos de 'close')."""
    v = points - pos
    d = np.linalg.norm(v, axis=1)
    cos = (v @ forward) / np.maximum(d, 1e-12)
    score = np.where((cos > np.cos(np.radians(cone_deg))) | (d < close), d, np.inf)
    k = int(np.argmin(score))
    return k if np.isfinite(score[k]) else None
