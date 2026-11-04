"""
La nave: una cámara libre en 6 ejes.

Velocidad logarítmica: la velocidad es proporcional a la distancia al objeto
más cercano. Lejos vas rápido, cerca vas despacio, y el viaje se siente igual
a cualquier escala (como en Space Engine o en "Powers of Ten").
"""
import numpy as np

from microscope.mathutil import orthonormalize, rotation


class Ship:
    def __init__(self, position, look_at):
        self.pos = np.array(position, dtype=float)
        self.vel = np.zeros(3)
        fwd = np.asarray(look_at, float) - self.pos
        fwd /= np.linalg.norm(fwd)
        up_hint = np.array([0.0, 0.0, 1.0]) if abs(fwd[2]) < 0.95 else np.array([0.0, 1.0, 0.0])
        right = np.cross(fwd, up_hint)
        right /= np.linalg.norm(right)
        up = np.cross(right, fwd)
        # columnas: derecha, arriba, atrás (convención OpenGL: se mira hacia -z)
        self.orient = np.column_stack([right, up, -fwd])
        self.speed_mult = 1.0
        self.mouse_sens = 0.0025

    @property
    def right(self):
        return self.orient[:, 0]

    @property
    def up(self):
        return self.orient[:, 1]

    @property
    def forward(self):
        return -self.orient[:, 2]

    def look(self, dx: float, dy: float):
        r = rotation(self.up, -dx * self.mouse_sens) @ rotation(self.right, dy * self.mouse_sens)
        self.orient = orthonormalize(r @ self.orient)

    def roll(self, angle: float):
        self.orient = orthonormalize(rotation(self.forward, angle) @ self.orient)

    def update(self, dt: float, move: np.ndarray, boost: bool, surface_dist: float,
               factor: float = 1.0):
        """move: (derecha, arriba, adelante) en [-1, 1]. factor: extra del piloto automático."""
        speed = (1.2 * max(surface_dist, 1e-9) * self.speed_mult * (5.0 if boost else 1.0)
                 * factor)
        target = (self.right * move[0] + self.up * move[1] + self.forward * move[2]) * speed
        # inercia suave
        k = 1 - np.exp(-dt * 6.0)
        self.vel += (target - self.vel) * k
        # nunca ir más rápido de lo que permite la distancia actual: así la nave
        # frena al acercarse aunque venga "lanzada" de un nivel superior
        v = np.linalg.norm(self.vel)
        if v > speed > 0:
            self.vel *= speed / v
        # nunca avanzar más de la mitad de la distancia en un frame (no "atravesar")
        step = self.vel * dt
        max_step = 0.5 * surface_dist
        n = np.linalg.norm(step)
        if n > max_step > 0:
            step *= max_step / n
        self.pos += step

    def rescale(self, center, scale):
        """Cambio de nivel: pasa a coordenadas del nuevo nivel."""
        self.pos = (self.pos - center) * scale
        self.vel = self.vel * scale

    def unscale(self, center, scale):
        self.pos = self.pos / scale + center
        self.vel = self.vel / scale
