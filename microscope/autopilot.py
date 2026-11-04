"""
Piloto automático.

  Enter        -> bajar un nivel (al objeto de la mira, o a uno elegido)
  Mayús+Enter  -> viaje completo hasta los quarks
  Retroceso    -> subir un nivel

El piloto no mueve la nave por su cuenta: decide hacia dónde mirar y le
dice al motor "avanza". La velocidad sigue siendo la logarítmica de la
nave, así que el viaje se siente igual que pilotando a mano.

Los textos (estado y avisos) salen del traductor del motor, así cambian con la tecla L.

Fases:
  turn   -> girar suavemente hasta mirar al destino
  fly    -> avanzar corrigiendo el rumbo (el destino puede moverse: quarks)
  face   -> tras subir de nivel, girarse para mirar el objeto del que venimos
"""
from typing import Optional

import numpy as np

from microscope.levels.base import Target
from microscope.mathutil import orthonormalize, rotation

TURN_RATE = 2.4          # rad/s al girar
ALIGNED_DEG = 1.5        # se considera alineado por debajo de este ángulo
SPEED_FACTOR = 2.0       # el piloto va el doble de rápido que la nave a mano


class Autopilot:
    def __init__(self, tr):
        self.tr = tr
        self.mode: Optional[str] = None      # "dive" | "full" | "ascend"
        self.phase: Optional[str] = None
        self.target: Optional[Target] = None
        self.point: Optional[np.ndarray] = None
        self.depth_at_start = 0
        self.exited_center = None
        self._clear_dist = 0.0
        self.message = ""                    # último aviso para el HUD (lo recoge main)

    # ------------------------------------------------------------ órdenes
    @property
    def active(self) -> bool:
        return self.mode is not None

    def dive(self, engine, full=False):
        self.mode = "full" if full else "dive"
        self._new_dive_target(engine)

    def ascend(self, engine):
        self.mode = "ascend"
        self.target = None
        self.depth_at_start = engine.level.depth
        p = engine.level.ascend_point(engine.ship.pos)
        self.point = p if p is not None else engine.home.copy()
        self.phase = "turn"

    def cancel(self):
        self.mode = self.phase = self.target = self.point = None

    @property
    def only_key(self):
        """El motor solo permite entrar en este hijo mientras el piloto vuela."""
        if self.mode in ("dive", "full") and self.target is not None:
            return self.target.key if self.target.key is not None else "__none__"
        if self.mode in ("ascend", "clear", "face"):
            return "__none__"            # subiendo: no entrar en nada por el camino
        return None

    def status(self) -> str:
        tr = self.tr
        if self.mode in ("dive", "full") and self.target:
            key = "microscope.autopilot.full" if self.mode == "full" else "microscope.autopilot.dive"
            return tr(key, target=self.target.label())
        if self.mode in ("ascend", "clear"):
            return tr("microscope.autopilot.ascend")
        return ""

    def flash_text(self) -> str:
        """Aviso breve al pulsar Enter: 'Piloto automático' o 'Viaje completo'."""
        return self.tr("microscope.autopilot.full_short" if self.mode == "full"
                       else "microscope.autopilot.dive_short")

    def target_point(self):
        if self.target is not None:
            return self.target.point()
        return self.point

    # ------------------------------------------------------------ lógica
    def _new_dive_target(self, engine):
        lvl = engine.level
        self.target = lvl.dive_target(engine.ship.pos, engine.ship.forward, engine.didactic)
        self.depth_at_start = lvl.depth
        self.phase = "turn"
        if self.target is None:
            self.cancel()

    def on_level_change(self, engine, entered: bool, exited_center=None):
        """Lo llama el motor al cambiar de nivel."""
        if not self.active:
            return
        if self.mode == "full" and entered:
            self._new_dive_target(engine)
            if self.target is not None:
                self.message = self.tr("microscope.autopilot.next", target=self.target.label())
        elif self.mode == "dive" and entered:
            self.cancel()
        elif self.mode == "ascend" and not entered:
            # ya fuera: seguir alejándose hasta no estar en la zona de entrada de
            # ningún otro objeto, y luego mirar hacia el objeto del que venimos
            self.mode, self.phase = "clear", "fly"
            away = engine.ship.pos - exited_center
            n = np.linalg.norm(away)
            away = away / n if n > 1e-12 else engine.ship.forward
            self.exited_center = np.array(exited_center, dtype=float)
            self.point = engine.ship.pos + away * max(n, 1e-12) * 50
            self._clear_dist = 1.6 * n
        else:
            self.cancel()

    def step(self, engine, dt):
        """Devuelve (mover, factor_velocidad). Gira la nave hacia el destino."""
        ship = engine.ship
        goal = self.target_point()
        if goal is None:
            self.cancel()
            return np.zeros(3), 1.0
        to_goal = goal - ship.pos
        dist = np.linalg.norm(to_goal)
        if dist < 1e-12:
            self.cancel()
            return np.zeros(3), 1.0
        angle = self._turn_towards(ship, to_goal / dist, dt)

        if self.phase == "turn":
            if angle < np.radians(ALIGNED_DEG):
                if self.mode == "face":
                    self.cancel()
                    return np.zeros(3), 1.0
                self.phase = "fly"
            return np.zeros(3), 1.0

        # fase "fly"
        if self.target is not None and self.target.key is None:
            # último nivel: llegar junto al quark y parar
            if dist < self.target.arrive_radius:
                self.message = self.tr("microscope.autopilot.end_of_journey")
                self.cancel()
                return np.zeros(3), 1.0
        if self.mode == "clear":
            far_enough = np.linalg.norm(ship.pos - self.exited_center) > self._clear_dist
            if engine.level.child_to_enter(ship.pos) is None and far_enough:
                self.mode, self.phase, self.point = "face", "turn", self.exited_center
                return np.zeros(3), 1.0
        if self.mode == "ascend" and engine.level.depth == 0:
            # en la raíz: volver a la vista panorámica, frenando al final
            if dist < 0.03 * np.linalg.norm(goal):
                self.mode, self.phase, self.point = "face", "turn", np.zeros(3)
                return np.zeros(3), 1.0
            brake = float(np.clip(dist / (0.4 * np.linalg.norm(goal)), 0.08, 1.0))
            return np.array([0.0, 0.0, brake]), SPEED_FACTOR
        # no avanzar mucho si aún no mira bien al destino (objetivo en movimiento)
        forward = 1.0 if angle < np.radians(20) else 0.2
        return np.array([0.0, 0.0, forward]), SPEED_FACTOR

    @staticmethod
    def _turn_towards(ship, direction, dt) -> float:
        f = ship.forward
        cosang = float(np.clip(f @ direction, -1, 1))
        angle = np.arccos(cosang)
        if angle < 1e-6:
            return 0.0
        axis = np.cross(f, direction)
        if np.linalg.norm(axis) < 1e-9:          # justo de espaldas: girar sobre "arriba"
            axis = ship.up
        step = min(angle, TURN_RATE * dt * (0.35 + 0.65 * min(angle / 0.5, 1.0)))
        ship.orient = orthonormalize(rotation(axis, step) @ ship.orient)
        return angle - step
