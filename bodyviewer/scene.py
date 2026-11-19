"""Las mallas del cuerpo listas para dibujar: un bloque por sistema, con la posición, la normal
y el número de pieza de cada vértice.

Coordenadas: BodyParts3D usa milímetros con z hacia arriba e y hacia atrás (la cara mira a -y).
Aquí se pasan a las de OpenGL (y hacia arriba, z hacia quien mira): (x, y, z) -> (x, z, -y), y se
centran en el cuerpo. Así la cámara por defecto, en +z, ve el cuerpo de frente.
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from dnainvent.anatomy.meshes import read_pack, vertex_normals
from dnainvent.anatomy.systems import SYSTEMS

COLORS = dict(SYSTEMS)


def to_gl(points_mm) -> np.ndarray:
    p = np.asarray(points_mm, np.float32)
    return np.stack([p[:, 0], p[:, 2], -p[:, 1]], axis=1)


def hex_rgb(color: str) -> tuple[float, float, float]:
    color = color.lstrip("#")
    return tuple(int(color[i:i + 2], 16) / 255 for i in (0, 2, 4))


@dataclass
class SystemMesh:
    system: str
    positions: np.ndarray          # float32 (N, 3), en mm, centrado
    normals: np.ndarray            # float32 (N, 3)
    piece: np.ndarray              # float32 (N,): número de pieza (índice en Scene.pieces)
    indices: np.ndarray            # uint32 (M,)


@dataclass
class Scene:
    meshes: list[SystemMesh]
    pieces: list[dict]             # {"id", "fma", "system"} por número de pieza
    center: np.ndarray             # el centro del cuerpo en coordenadas de BodyParts3D (gl)
    radius: float
    lo: np.ndarray                 # caja de cada pieza (P, 3), ya centrada
    hi: np.ndarray
    by_id: dict = field(default_factory=dict)          # id de pieza -> número

    @property
    def systems(self) -> list[str]:
        return [m.system for m in self.meshes]

    @property
    def triangles(self) -> int:
        return sum(len(m.indices) // 3 for m in self.meshes)

    def bounds(self, numbers) -> tuple[np.ndarray, np.ndarray] | None:
        numbers = list(numbers)
        if not numbers:
            return None
        return self.lo[numbers].min(axis=0), self.hi[numbers].max(axis=0)


def load_scene(meshes_dir) -> Scene:
    """Lee meshes/<sistema>.dnab en el orden de las capas (systems.SYSTEMS)."""
    meshes_dir = Path(meshes_dir)
    packs = [(key, meshes_dir / f"{key}.dnab") for key, _ in SYSTEMS]
    packs = [(key, read_pack(path)) for key, path in packs if path.is_file()]
    if not packs:
        raise FileNotFoundError(f"No hay mallas en {meshes_dir}")
    pieces, by_id, lo, hi = [], {}, [], []
    raw = []
    for key, pack in packs:
        positions, normals, numbers, indices, offset = [], [], [], [], 0
        for i, entry in enumerate(pack.pieces):
            vertices, faces = pack.mesh(i)
            vertices = to_gl(vertices)
            number = by_id.get(entry["id"])
            if number is None:
                number = by_id[entry["id"]] = len(pieces)
                pieces.append({"id": entry["id"], "fma": entry["fma"], "system": key})
                lo.append(vertices.min(axis=0))
                hi.append(vertices.max(axis=0))
            else:                     # una pieza partida en dos trozos (más de 65.535 vértices)
                lo[number] = np.minimum(lo[number], vertices.min(axis=0))
                hi[number] = np.maximum(hi[number], vertices.max(axis=0))
            positions.append(vertices)
            normals.append(vertex_normals(vertices, faces))
            numbers.append(np.full(len(vertices), number, np.float32))
            indices.append(faces.reshape(-1).astype(np.uint32) + offset)
            offset += len(vertices)
        raw.append((key, np.concatenate(positions), np.concatenate(normals), np.concatenate(numbers),
                    np.concatenate(indices)))
    lo, hi = np.array(lo, np.float32), np.array(hi, np.float32)
    center = (lo.min(axis=0) + hi.max(axis=0)) / 2
    radius = float(np.linalg.norm(hi.max(axis=0) - lo.min(axis=0)) / 2)
    meshes = [SystemMesh(key, (p - center).astype(np.float32), n.astype(np.float32), num, idx)
              for key, p, n, num, idx in raw]
    return Scene(meshes, pieces, center, radius, lo - center, hi - center, by_id)
