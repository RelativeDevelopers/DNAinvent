"""Estructura de proteína (PDB) -> molécula 3D para el microscopio, con todos sus átomos.

Tiene la misma forma que Molecule (scene/molecule.py): los mismos campos con los mismos
nombres. Así los niveles del microscopio que reciben una molécula (átomo, núcleo, nucleón)
sirven sin cambios para el ADN y para las proteínas. Las equivalencias:

    Molecule (ADN)                  ProteinMolecule
    sequence   bases del tramo      aminoácidos del tramo (1 letra)
    pair_index nº de par            nº de aminoácido en la cadena completa (desde 0)
    resn       base ("A")           aminoácido en 3 letras ("MET")
    part       base/sugar/phosph.   backbone (esqueleto) / sidechain (cadena lateral)
    strand     hebra 1 o 2          siempre 1 (una sola cadena)
    hbonds     puentes A-T, G-C     vacío por ahora

Los hidrógenos se calculan sobre la cadena completa (el H del enlace peptídico necesita al
aminoácido anterior) y después se recorta el tramo pedido, como en el ADN.

Las posiciones van en nm y centradas en el origen, como en Molecule. Los enlaces salen de
la tabla de cada aminoácido (no de distancias), más los peptídicos C(i)-N(i+1) y los
puentes disulfuro S-S.
"""
from dataclasses import dataclass, field
from typing import ClassVar

import numpy as np

from ..errors import InvalidSequenceError
from ..protein.geometry import distance
from ..protein.hydrogens import DISULFIDE_MAX, add_hydrogens, templates
from ..protein.pdb import Residue, Structure
from .limits import window

ANGSTROM_TO_NM = 0.1
BACKBONE = frozenset({"N", "CA", "C", "O", "OXT", "H", "H1", "H2", "H3", "HA", "HA2", "HA3"})
PEPTIDE_BOND_MAX = 2.0     # Å


@dataclass
class ProteinMolecule:
    kind: ClassVar[str] = "protein"        # Molecule (ADN) no lo tiene: se toma "dna"

    sequence: str                    # aminoácidos del tramo (1 letra)
    start: int                       # índice del primer aminoácido del tramo en la cadena
    length: int                      # aminoácidos de la cadena completa
    positions: np.ndarray            # (N, 3) en nm, centrada en el origen
    elements: list                   # símbolo por átomo
    names: list                      # nombre PDB del átomo ("CA", "HB2"...)
    strand: np.ndarray               # 1 (compatibilidad con Molecule)
    pair_index: np.ndarray           # nº de aminoácido en la cadena completa (desde 0)
    resn: list                       # aminoácido en 3 letras
    part: list                       # "backbone" | "sidechain"
    bonds: np.ndarray                # (M, 2) enlaces covalentes
    hbonds: np.ndarray               # (K, 2) vacío por ahora
    calculated_h: np.ndarray         # True si el H es calculado (no venía en el archivo)
    confidence: np.ndarray | None = None      # pLDDT por átomo (solo en predicciones)
    residue_numbers: list = field(default_factory=list)   # numeración del archivo, por residuo
    axis: np.ndarray = field(default_factory=lambda: np.array([0.0, 0.0, 1.0]))

    @property
    def n_atoms(self) -> int:
        return len(self.elements)

    @property
    def end(self) -> int:
        return self.start + len(self.sequence)

    @property
    def residue_index(self) -> np.ndarray:
        """Nombre más claro de pair_index para las proteínas."""
        return self.pair_index

    def radius_vdw(self) -> np.ndarray:
        from .elements import atom_type      # necesita AtomicWorld: solo se importa si se usa
        return np.array([atom_type(e).vdw_nm for e in self.elements], dtype=np.float32)

    def bonded_to(self, i: int) -> list[int]:
        rows = np.where((self.bonds[:, 0] == i) | (self.bonds[:, 1] == i))[0]
        return [int(b if a == i else a) for a, b in self.bonds[rows]]


def build_protein(structure: Structure, start: int = 0, count: int | None = None,
                  chain: str | None = None, hydrogens: bool = True,
                  origin: tuple[float, float, float] | None = None) -> ProteinMolecule:
    """Todos los átomos de los aminoácidos [start, start + count) de una cadena
    (por defecto, la primera), con los hidrógenos calculados si hydrogens=True.

    origin (Å): punto que queda en el centro. Por defecto, el centro del propio tramo. Para
    comparar dos proteínas ya superpuestas se les da el mismo origen, así comparten
    coordenadas (y la nave ve las dos en el mismo sitio)."""
    residues = structure.residues(chain)
    if not residues:
        raise InvalidSequenceError("empty_sequence", kind="protein")
    original = {(r.resi, r.icode): set(r.atoms) for r in residues}
    if hydrogens:
        residues = add_hydrogens(residues).residues
    start, end = window(len(residues), start, count)
    chosen = residues[start:end]

    pos, els, names, pidx, resn, part, calc, conf = [], [], [], [], [], [], [], []
    index = {}                                  # (k, nombre) -> nº de átomo
    for k, residue in enumerate(chosen, start=start):
        known = original[(residue.resi, residue.icode)]
        for name, atom in residue.atoms.items():
            index[(k, name)] = len(els)
            pos.append(atom.xyz)
            els.append(atom.element)
            names.append(name)
            pidx.append(k)
            resn.append(residue.resn)
            part.append("backbone" if name in BACKBONE else "sidechain")
            calc.append(atom.element == "H" and name not in known)
            conf.append(atom.bfactor)

    positions = np.array(pos, dtype=np.float64) * ANGSTROM_TO_NM
    positions -= (positions.mean(axis=0) if origin is None
                  else np.array(origin, dtype=np.float64) * ANGSTROM_TO_NM)
    bonds = _bonds(chosen, start, index)
    return ProteinMolecule(
        sequence="".join(r.letter for r in chosen), start=start, length=len(residues),
        positions=positions.astype(np.float32), elements=els, names=names,
        strand=np.ones(len(els), dtype=np.int8), pair_index=np.array(pidx), resn=resn,
        part=part, bonds=np.array(bonds, dtype=np.int32).reshape(-1, 2),
        hbonds=np.zeros((0, 2), dtype=np.int32), calculated_h=np.array(calc, dtype=bool),
        confidence=np.array(conf, dtype=np.float32) if structure.is_prediction else None,
        residue_numbers=[r.resi for r in chosen],
    )


def _bonds(residues: list[Residue], start: int, index: dict) -> list[list[int]]:
    data = templates()["residues"]
    bonds = []

    def link(k1, n1, k2, n2):
        a, b = index.get((k1, n1)), index.get((k2, n2))
        if a is not None and b is not None:
            bonds.append([a, b])

    for k, residue in enumerate(residues, start=start):
        template = data.get(residue.resn)
        if template:
            for n1, n2 in template["bonds"]:
                link(k, n1, k, n2)
        link(k, "C", k, "OXT")
        for h in ("H1", "H2", "H3"):
            link(k, "N", k, h)
        if k > start:
            previous = residues[k - start - 1]
            c, n = previous.atoms.get("C"), residue.atoms.get("N")
            if c is not None and n is not None and distance(c.xyz, n.xyz) <= PEPTIDE_BOND_MAX:
                link(k - 1, "C", k, "N")

    sulfurs = [(k, r.atoms["SG"]) for k, r in enumerate(residues, start=start)
               if r.resn == "CYS" and "SG" in r.atoms]
    for i, (k1, s1) in enumerate(sulfurs):
        for k2, s2 in sulfurs[i + 1:]:
            if distance(s1.xyz, s2.xyz) <= DISULFIDE_MAX:
                link(k1, "SG", k2, "SG")
    return bonds
