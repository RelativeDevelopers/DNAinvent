"""Añadir los hidrógenos que faltan en una estructura de proteína.

AlphaFold DB y ESMFold solo dan los átomos pesados (C, N, O, S). Los hidrógenos se calculan
con la tabla templates/amino_acids.json (generada con PyMOL: ver
tools/generate_amino_acid_templates.py): para cada hidrógeno, a qué átomo va unido y sus
coordenadas internas respecto a dos átomos pesados vecinos (ver geometry.place).

Simplificaciones (las mismas que en la tabla):
  - pH 7: Lys y Arg con carga +, Asp y Glu con carga - (sin H ácido).
  - Histidina con el H en Nε (HE2), una de sus dos formas neutras.
  - Grupos que giran (-CH3, -OH, -NH3+, -SH) en posición alternada, sin optimizar puentes
    de hidrógeno con los vecinos (lo que hacen programas como Reduce). Solo una corrección:
    si el H de un -OH o -SH choca con otro átomo (< 1,6 Å), se gira a la otra posición
    alternada más despejada (hay tres, separadas 120°).
  - Extremo N cargado (-NH3+; en la prolina, -NH2+). El extremo C (COO-) no lleva H.
  - Una cisteína unida a otra (puente disulfuro S-S) no lleva el H del azufre.

Son hidrógenos calculados por geometría, no medidos: los de posición fija quedan a
~0,1 Å de los reales; los de grupos que giran pueden estar girados.
"""

import json
from dataclasses import dataclass, field, replace
from functools import lru_cache
from pathlib import Path

from .geometry import distance, place
from .pdb import Atom, Residue

TEMPLATE_FILE = Path(__file__).resolve().parent / "templates" / "amino_acids.json"
PEPTIDE_BOND_MAX = 2.0         # Å: C(i)-N(i+1) más largo que esto = la cadena está cortada
DISULFIDE_MAX = 2.5            # Å: SG-SG más cerca que esto = puente disulfuro
CLASH = 1.6                    # Å: un H más cerca que esto de otro átomo no enlazado choca


@lru_cache(maxsize=1)
def templates() -> dict:
    return json.loads(TEMPLATE_FILE.read_text(encoding="utf-8"))


@dataclass
class HydrogenResult:
    residues: list[Residue]
    added: int = 0
    skipped: list[tuple[int, str]] = field(default_factory=list)   # (residuo, nombre del H)


def add_hydrogens(residues: list[Residue], keep_existing: bool = True) -> HydrogenResult:
    """Devuelve una copia de los residuos con sus hidrógenos.

    keep_existing: si la estructura ya trae hidrógenos (p. ej. una de RMN), se conservan y
    solo se añaden los que falten. Con False se quitan todos y se recalculan.

    Si faltan los átomos de referencia de algún H (estructuras experimentales incompletas),
    ese H se salta y se anota en .skipped.
    """
    data = templates()
    disulfide = _disulfide_cysteines(residues)
    result = HydrogenResult([])
    rotors = []                               # (residuo, spec, coords) de cada -OH / -SH

    for i, residue in enumerate(residues):
        atoms = {name: atom for name, atom in residue.atoms.items()
                 if keep_existing or atom.element != "H"}
        new = Residue(residue.chain, residue.resi, residue.icode, residue.resn, dict(atoms))
        result.residues.append(new)
        template = data["residues"].get(residue.resn)
        if template is None:                  # ligando o aminoácido no estándar
            continue

        previous = residues[i - 1] if i > 0 else None
        linked = previous is not None and _peptide_bonded(previous, residue)
        coords = {name: atom.xyz for name, atom in atoms.items()}
        if linked:
            coords["-C"] = previous.atoms["C"].xyz

        specs = list(template["hydrogens"])
        if not linked:                        # extremo N: -NH3+ en lugar del H del enlace
            specs = [h for h in specs if h["name"] != "H"]
            kind = "proline" if residue.resn == "PRO" else "regular"
            specs += data["n_terminus"][kind]
        if (residue.chain, residue.resi, residue.icode) in disulfide:
            specs = [h for h in specs if h["name"] != "HG"]

        for h in specs:
            if h["name"] in atoms:
                continue
            r2, r3 = h["refs"]
            try:
                position = place(coords[r3], coords[r2], coords[h["parent"]],
                                 h["r"], h["angle"], h["dihedral"])
            except KeyError:
                result.skipped.append((residue.resi, h["name"]))
                continue
            parent = atoms[h["parent"]]
            new.atoms[h["name"]] = replace(parent, name=h["name"], element="H",
                                           x=position[0], y=position[1], z=position[2])
            result.added += 1
            if h["rotates"] and parent.element in ("O", "S"):
                rotors.append((new, h, coords))

    _avoid_clashes(rotors, result.residues)
    return result


def _avoid_clashes(rotors, residues: list[Residue]) -> None:
    """Gira el H de un -OH o -SH a otra posición alternada si choca con algún átomo."""
    if not rotors:
        return
    everything = [(r, name, atom.xyz) for r in residues for name, atom in r.atoms.items()]

    def clearance(residue, h, position):
        ignore = {h["name"], h["parent"], h["refs"][0]}
        return min((distance(position, xyz) for r, name, xyz in everything
                    if not (r is residue and name in ignore)), default=float("inf"))

    for residue, h, coords in rotors:
        atom = residue.atoms[h["name"]]
        if clearance(residue, h, atom.xyz) >= CLASH:
            continue
        r2, r3 = h["refs"]
        candidates = [place(coords[r3], coords[r2], coords[h["parent"]], h["r"], h["angle"],
                            h["dihedral"] + turn) for turn in (0.0, 120.0, -120.0)]
        best = max(candidates, key=lambda p: clearance(residue, h, p))
        residue.atoms[h["name"]] = replace(atom, x=best[0], y=best[1], z=best[2])
        everything = [(r, n, a.xyz) for r in residues for n, a in r.atoms.items()]


def is_rotatable(resn: str, name: str) -> bool:
    """True si el H pertenece a un grupo que gira (su posición exacta es una elección)."""
    data = templates()
    for h in data["residues"].get(resn, {}).get("hydrogens", []):
        if h["name"] == name:
            return h["rotates"]
    return name in ("H1", "H2", "H3")         # extremo N


def _peptide_bonded(previous: Residue, residue: Residue) -> bool:
    c, n = previous.atoms.get("C"), residue.atoms.get("N")
    return c is not None and n is not None and distance(c.xyz, n.xyz) <= PEPTIDE_BOND_MAX


def _disulfide_cysteines(residues: list[Residue]) -> set:
    sulfurs = [(r, r.atoms["SG"]) for r in residues if r.resn == "CYS" and "SG" in r.atoms]
    bonded = set()
    for k, (r1, s1) in enumerate(sulfurs):
        for r2, s2 in sulfurs[k + 1:]:
            if distance(s1.xyz, s2.xyz) <= DISULFIDE_MAX:
                bonded |= {(r1.chain, r1.resi, r1.icode), (r2.chain, r2.resi, r2.icode)}
    return bonded
