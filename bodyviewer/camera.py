"""La cámara del visor: gira alrededor de un punto (el cuerpo o la pieza elegida)."""

import math

import numpy as np


def perspective(fov_y: float, aspect: float, near: float, far: float) -> np.ndarray:
    f = 1 / math.tan(fov_y / 2)
    m = np.zeros((4, 4), np.float32)
    m[0, 0], m[1, 1] = f / aspect, f
    m[2, 2], m[2, 3] = (far + near) / (near - far), 2 * far * near / (near - far)
    m[3, 2] = -1
    return m


def look_at(eye, target, up=(0, 1, 0)) -> np.ndarray:
    eye, target, up = (np.asarray(v, np.float64) for v in (eye, target, up))
    f = target - eye
    f /= np.linalg.norm(f)
    s = np.cross(f, up)
    s /= np.linalg.norm(s)
    u = np.cross(s, f)
    m = np.identity(4)
    m[0, :3], m[1, :3], m[2, :3] = s, u, -f
    m[:3, 3] = -m[:3, :3] @ eye
    return m.astype(np.float32)


class OrbitCamera:
    FOV = math.radians(40)
    MIN_PITCH, MAX_PITCH = -1.45, 1.45

    def __init__(self, radius: float):
        self.home_radius = radius
        self.reset()

    def reset(self):
        self.target = np.zeros(3, np.float32)
        self.yaw, self.pitch = 0.0, 0.05
        self.distance = self.fit_distance(self.home_radius)

    def fit_distance(self, radius: float) -> float:
        return max(radius, 5.0) / math.sin(self.FOV / 2) * 1.05

    def frame(self, lo, hi):
        """Que la caja lo-hi quepa en la pantalla, mirando a su centro."""
        lo, hi = np.asarray(lo, np.float32), np.asarray(hi, np.float32)
        self.target = (lo + hi) / 2
        self.distance = self.fit_distance(float(np.linalg.norm(hi - lo)) / 2)

    @property
    def eye(self) -> np.ndarray:
        cp = math.cos(self.pitch)
        offset = np.array([math.sin(self.yaw) * cp, math.sin(self.pitch), math.cos(self.yaw) * cp])
        return self.target + offset * self.distance

    def orbit(self, dx: float, dy: float):
        self.yaw -= dx * 0.008
        self.pitch = float(np.clip(self.pitch - dy * 0.008, self.MIN_PITCH, self.MAX_PITCH))

    def zoom(self, steps: float):
        self.distance = float(np.clip(self.distance * (0.88 ** steps), 30.0, self.home_radius * 8))

    def pan(self, dx: float, dy: float, height_px: int):
        """Mover el punto al que se mira, en el plano de la pantalla (px -> mm)."""
        view = self.view()
        right, up = view[0, :3], view[1, :3]
        mm_per_px = 2 * self.distance * math.tan(self.FOV / 2) / max(height_px, 1)
        self.target = (self.target - right * dx * mm_per_px - up * dy * mm_per_px).astype(np.float32)

    def view(self) -> np.ndarray:
        return look_at(self.eye, self.target)

    def projection(self, aspect: float) -> np.ndarray:
        near = max(self.distance * 0.01, 1.0)
        return perspective(self.FOV, aspect, near, self.distance + self.home_radius * 4)
