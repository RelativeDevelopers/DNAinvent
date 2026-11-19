"""Las mallas del cuerpo, empaquetadas: un archivo por sistema (huesos, músculos...), que leen
igual la consola (visor de escritorio) y la web (WebGL).

Formato «DNAB» versión 1, todo en little-endian:

    b"DNAB"  uint32 versión  uint32 bytes de la cabecera  cabecera JSON (UTF-8, rellena con
    espacios hasta múltiplo de 4)  posiciones  índices

Cabecera:
    {"system": "skeletal", "origin": [x, y, z], "step": mm por unidad,
     "pieces": [{"id": "FJ3155", "fma": "FMA24474", "vertex": [primero, cuántos],
                 "index": [primero, cuántos]}, ...]}

Posiciones: uint16 x, y, z por vértice (mm = origin + valor × step; el mismo origen y paso en
todos los archivos, así las piezas encajan). Índices: uint16, contados desde el primer vértice
de su pieza (ninguna pieza pasa de 65.535 vértices; si pasara, se parte). Las normales no se
guardan: el visor las calcula (vertex_normals). Coordenadas de BodyParts3D: milímetros, z hacia
arriba (la cabeza), y hacia atrás.
"""

import json
import struct
from dataclasses import dataclass
from pathlib import Path

import numpy as np

MAGIC = b"DNAB"
VERSION = 1
MAX_VERTICES = 65535


def simplify(vertices, faces, ratio: float, min_faces: int = 300):
    """Menos triángulos, la misma forma (cuádricas de error, con fast-simplification). Las
    piezas pequeñas se quedan como están. Sin la librería, no se simplifica."""
    target = max(min_faces, int(len(faces) * ratio))
    if ratio >= 1 or target >= len(faces):
        return vertices, faces
    try:
        import fast_simplification
    except ImportError:
        return vertices, faces
    v, f = fast_simplification.simplify(np.ascontiguousarray(vertices, np.float32),
                                        np.ascontiguousarray(faces, np.int32), target_count=target)
    return np.asarray(v, np.float32), np.asarray(f, np.int32)


def quantize(vertices, origin, step):
    q = np.round((np.asarray(vertices, np.float64) - origin) / step)
    return np.clip(q, 0, 65535).astype(np.uint16)


def weld(quantized, faces):
    """Junta los vértices que caen en el mismo sitio (los OBJ repiten muchos) y quita los
    triángulos que se quedan sin área."""
    unique, inverse = np.unique(quantized, axis=0, return_inverse=True)
    f = inverse.reshape(-1)[faces]
    keep = (f[:, 0] != f[:, 1]) & (f[:, 1] != f[:, 2]) & (f[:, 0] != f[:, 2])
    return unique, f[keep]


def split(quantized, faces):
    """Trozos de como mucho MAX_VERTICES vértices (para que los índices quepan en uint16)."""
    if len(quantized) <= MAX_VERTICES:
        return [(quantized, faces)]
    out, start = [], 0
    while start < len(faces):
        used, chunk = {}, []
        for tri in faces[start:]:
            new = [v for v in tri if v not in used]
            if len(used) + len(new) > MAX_VERTICES:
                break
            for v in new:
                used[v] = len(used)
            chunk.append([used[v] for v in tri])
            start += 1
        order = np.array(sorted(used, key=used.get))
        out.append((quantized[order], np.array(chunk, np.int64)))
    return out


def write_pack(path, system: str, pieces, origin, step) -> dict:
    """pieces: [(id, fma, posiciones uint16 N×3, triángulos M×3)]. Devuelve la cabecera."""
    header = {"system": system, "origin": [float(x) for x in origin], "step": float(step), "pieces": []}
    positions, indices, nv, ni = [], [], 0, 0
    for piece_id, fma, q, faces in pieces:
        for part_q, part_f in split(q, faces):
            header["pieces"].append({"id": piece_id, "fma": fma, "vertex": [nv, len(part_q)],
                                     "index": [ni, part_f.size]})
            positions.append(np.ascontiguousarray(part_q, np.uint16))
            indices.append(np.ascontiguousarray(part_f, np.uint16).reshape(-1))
            nv += len(part_q)
            ni += part_f.size
    blob = json.dumps(header, separators=(",", ":")).encode("utf-8")
    blob += b" " * (-len(blob) % 4)
    with open(path, "wb") as out:
        out.write(MAGIC + struct.pack("<II", VERSION, len(blob)) + blob)
        for p in positions:
            out.write(p.tobytes())
        for i in indices:
            out.write(i.tobytes())
        if ni % 2:
            out.write(b"\0\0")                     # que el archivo acabe alineado a 4
    return header


@dataclass
class Pack:
    header: dict
    positions: np.ndarray          # uint16, (N, 3)
    indices: np.ndarray            # uint16, (M,)

    @property
    def pieces(self) -> list[dict]:
        return self.header["pieces"]

    def mesh(self, i: int):
        """Vértices (mm, float32) y triángulos (desde 0) de la entrada i."""
        piece = self.pieces[i]
        v0, nv = piece["vertex"]
        i0, ni = piece["index"]
        origin = np.array(self.header["origin"], np.float32)
        vertices = origin + self.positions[v0:v0 + nv].astype(np.float32) * self.header["step"]
        return vertices, self.indices[i0:i0 + ni].astype(np.int32).reshape(-1, 3)


def read_pack(path) -> Pack:
    data = Path(path).read_bytes()
    if data[:4] != MAGIC:
        raise ValueError(f"{path} no es un paquete de mallas de DNAinvent")
    version, size = struct.unpack_from("<II", data, 4)
    if version != VERSION:
        raise ValueError(f"{path}: versión {version} (se esperaba {VERSION})")
    header = json.loads(data[12:12 + size].decode("utf-8"))
    nv = sum(p["vertex"][1] for p in header["pieces"])
    ni = sum(p["index"][1] for p in header["pieces"])
    start = 12 + size
    positions = np.frombuffer(data, np.uint16, nv * 3, start).reshape(-1, 3)
    indices = np.frombuffer(data, np.uint16, ni, start + nv * 6)
    return Pack(header, positions, indices)


def vertex_normals(vertices, faces):
    """Normales suaves: la media de las de los triángulos que tocan cada vértice."""
    v = np.asarray(vertices, np.float32)
    f = np.asarray(faces)
    n = np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]])
    out = np.zeros_like(v)
    for k in range(3):
        np.add.at(out, f[:, k], n)
    length = np.linalg.norm(out, axis=1, keepdims=True)
    return out / np.where(length == 0, 1, length)
