"""Matrices de cámara (sin dependencias aparte de NumPy)."""
import numpy as np


def perspective(fov_y_deg: float, aspect: float, near: float, far: float) -> np.ndarray:
    f = 1.0 / np.tan(np.radians(fov_y_deg) / 2)
    m = np.zeros((4, 4))
    m[0, 0] = f / aspect
    m[1, 1] = f
    m[2, 2] = (far + near) / (near - far)
    m[2, 3] = 2 * far * near / (near - far)
    m[3, 2] = -1
    return m


def view_matrix(position: np.ndarray, orientation: np.ndarray) -> np.ndarray:
    """orientation: columnas = (derecha, arriba, atrás) en coordenadas del mundo."""
    m = np.eye(4)
    m[:3, :3] = orientation.T
    m[:3, 3] = -orientation.T @ position
    return m


def rotation(axis: np.ndarray, angle: float) -> np.ndarray:
    axis = axis / np.linalg.norm(axis)
    x, y, z = axis
    c, s = np.cos(angle), np.sin(angle)
    C = 1 - c
    return np.array([[c + x * x * C, x * y * C - z * s, x * z * C + y * s],
                     [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
                     [z * x * C - y * s, z * y * C + x * s, c + z * z * C]])


def orthonormalize(m: np.ndarray) -> np.ndarray:
    u, _, vt = np.linalg.svd(m)
    return u @ vt


def gl(m: np.ndarray) -> bytes:
    """Matriz NumPy (convención fila-columna) -> bytes en el orden que espera GLSL."""
    return m.T.astype("f4").tobytes()
