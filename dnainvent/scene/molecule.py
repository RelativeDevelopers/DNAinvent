"""Secuencia "ATGC..." -> molécula de ADN 3D (forma B) con todos sus átomos.

Cada par de bases se toma de una plantilla atómica real (templates/b_dna.json, generada con
PyMOL: tools/generate_templates.py) y se coloca aplicando k veces el paso helicoidal.

El paso es el mismo que usa dnainvent/helix3d.py (y por tanto la web 2D): un giro de
360° / 10,5 ≈ 34,3° y una subida de 3,38 Å por par, alrededor del eje de la plantilla.
Así la hélice atómica y la del modelo sencillo tienen la misma forma.

Las secuencias largas no se construyen enteras: se elige un tramo (start, count). Los
números de par (pair_index) son los de la secuencia completa, para poder señalarlos.
"""
import json
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import numpy as np

from ..constants import COMPLEMENT
from ..helix3d import BASE_PAIRS_PER_TURN, RISE
from ..validation import validate_dna
from .elements import atom_type
from .limits import RECOMMENDED_PAIRS, window  # noqa: F401  (se reexportan)

TEMPLATE_FILE = Path(__file__).resolve().parent / "templates" / "b_dna.json"
ANGSTROM_TO_NM = 0.1
HBOND_MIN_NM, HBOND_MAX_NM = 0.26, 0.33     # N/O···N/O de un puente de hidrógeno


@dataclass
class Molecule:
    sequence: str                    # bases del tramo construido (hebra 1, 5'->3')
    start: int                       # posición del primer par en la secuencia completa
    length: int                      # longitud de la secuencia completa
    positions: np.ndarray            # (N, 3) en nm, centrada en el origen
    elements: list                   # símbolo por átomo
    names: list                      # nombre PDB del átomo ("N7", "C1'"...)
    strand: np.ndarray               # 1 o 2
    pair_index: np.ndarray           # número de par en la secuencia completa (desde 0)
    resn: list                       # base del nucleótido al que pertenece el átomo
    part: list                       # "base" | "sugar" | "phosphate"
    bonds: np.ndarray                # (M, 2) enlaces covalentes
    hbonds: np.ndarray               # (K, 2) puentes de hidrógeno entre bases
    axis: np.ndarray = field(default_factory=lambda: np.array([0.0, 0.0, 1.0]))

    @property
    def n_atoms(self) -> int:
        return len(self.elements)

    @property
    def end(self) -> int:
        return self.start + len(self.sequence)

    def radius_vdw(self) -> np.ndarray:
        return np.array([atom_type(e).vdw_nm for e in self.elements], dtype=np.float32)

    def bonded_to(self, i: int) -> list[int]:
        rows = np.where((self.bonds[:, 0] == i) | (self.bonds[:, 1] == i))[0]
        return [int(b if a == i else a) for a, b in self.bonds[rows]]


@lru_cache(maxsize=1)
def _templates():
    data = json.loads(TEMPLATE_FILE.read_text(encoding="utf-8"))
    return data["pairs"], np.array(data["step"]["R"]), np.array(data["step"]["t"])


def helical_step() -> tuple[np.ndarray, np.ndarray]:
    """(R, t) del paso entre un par y el siguiente: x_k+1 = R · x_k + t.

    La plantilla trae el paso de PyMOL (10 pares por vuelta). Se conserva su eje (el
    punto c por el que pasa, que no es el origen) y su sentido, y se cambian el giro
    y la subida por los de DNAinvent: es un movimiento de tornillo alrededor de ese eje.
    """
    _, R0, t0 = _templates()
    c = np.linalg.solve((np.eye(3) - R0)[:2, :2], t0[:2])       # punto del eje
    sign = np.sign(R0[1, 0])                                     # sentido del giro
    theta = sign * np.radians(360.0 / BASE_PAIRS_PER_TURN)
    R = np.array([[np.cos(theta), -np.sin(theta), 0.0],
                  [np.sin(theta), np.cos(theta), 0.0],
                  [0.0, 0.0, 1.0]])
    t = np.zeros(3)
    t[:2] = (np.eye(2) - R[:2, :2]) @ c
    t[2] = np.sign(t0[2]) * RISE
    return R, t


def build_molecule(dna: str, start: int = 0, count: int | None = None) -> Molecule:
    """Todos los átomos del tramo [start, start + count) de la doble hélice de 'dna'."""
    dna = validate_dna(dna)
    start, end = window(len(dna), start, count)
    pairs, _, _ = _templates()
    R, t = helical_step()

    pos, els, names, strand, pidx, resn, part = [], [], [], [], [], [], []
    bonds = []
    o3_prev = None          # O3' de la hebra 1 en el par anterior
    p_prev = None           # P de la hebra 2 en el par anterior

    Rk, tk = np.eye(3), np.zeros(3)
    for k in range(start, end):
        base = dna[k]
        tpl = pairs[base + COMPLEMENT[base]]
        offset = len(els)
        local = {}
        for i, a in enumerate(tpl["atoms"]):
            pos.append(Rk @ np.array(a["xyz"]) + tk)
            els.append(a["el"])
            names.append(a["name"])
            strand.append(a["strand"])
            pidx.append(k)
            resn.append(a["resn"])
            part.append(a["part"])
            local[(a["strand"], a["name"])] = offset + i
        bonds.extend([offset + i, offset + j] for i, j in tpl["bonds"])
        # esqueleto entre pares: hebra 1 (5'->3' en el sentido de k) O3'(k-1) — P(k);
        # hebra 2 (antiparalela) O3'(k) — P(k-1)
        if o3_prev is not None:
            bonds.append([o3_prev, local[(1, "P")]])
            bonds.append([local[(2, "O3'")], p_prev])
        o3_prev, p_prev = local[(1, "O3'")], local[(2, "P")]
        Rk, tk = R @ Rk, R @ tk + t

    positions = np.array(pos) * ANGSTROM_TO_NM
    positions -= positions.mean(axis=0)
    strand, pidx = np.array(strand), np.array(pidx)
    hbonds = find_hbonds(positions, els, strand, pidx, part)
    return Molecule(dna[start:end], start, len(dna), positions.astype(np.float32), els, names,
                    strand, pidx, resn, part, np.array(bonds, dtype=np.int32).reshape(-1, 2),
                    np.array(hbonds, dtype=np.int32).reshape(-1, 2))


def find_hbonds(positions, els, strand, pidx, part) -> list[list[int]]:
    """Puentes de hidrógeno N/O···N/O entre las dos bases de cada par (2,6–3,3 Å).

    No se escriben a mano: salen de la geometría. A–T da 2 y G–C da 3.
    """
    out = []
    for k in np.unique(pidx):
        members = np.where(pidx == k)[0]
        polar = [i for i in members if part[i] == "base" and els[i] in ("N", "O")]
        s1 = [i for i in polar if strand[i] == 1]
        s2 = [i for i in polar if strand[i] == 2]
        for i in s1:
            for j in s2:
                if HBOND_MIN_NM <= np.linalg.norm(positions[i] - positions[j]) <= HBOND_MAX_NM:
                    out.append([int(i), int(j)])
    return out
