"""Moléculas dentro de una proteína (fármacos, cofactores...) y cómo encajan en ella.

En una estructura experimental del PDB, todo lo que no es la cadena de aminoácidos va en
líneas HETATM: el fármaco, pero también el agua, el hemo de la hemoglobina, iones y restos
del cristal (sulfato, glicerol...). Aquí se separan:

    drug        lo interesante: una molécula pequeña que no es de las de siempre
    cofactor    moléculas que la proteína necesita para funcionar (hemo, NAD, ATP...) o gases
                unidos a ellas (CO, O2)
    additive    iones y restos del experimento (sulfato, glicerol, sodio...)
    (el agua se descarta)

Contactos: para cada átomo del fármaco, los átomos de la proteína a menos de 4 Å. No es una
simulación: son las distancias medidas en el cristal, y el tipo se deduce de los elementos y
la distancia, como hacen los visores de estructuras:

    covalent      < 2,0 Å entre átomos pesados: el fármaco está unido químicamente (el
                  voxelotor forma un enlace con el primer aminoácido de la cadena alfa)
    hbond         N u O con N u O a <= 3,5 Å: puente de hidrógeno
    hydrophobic   C o S con C o S a <= 4,0 Å: zonas «aceitosas» que se juntan
    contact       cualquier otro par a <= 4,0 Å

Los hidrógenos no cuentan (las estructuras de cristal casi nunca los traen).
"""

import math
from dataclasses import dataclass

from .pdb import Atom, Structure

CONTACT_CUTOFF = 4.0
HBOND_CUTOFF = 3.5
COVALENT_CUTOFF = 2.0
MIN_DRUG_ATOMS = 6            # menos átomos pesados: ion o resto pequeño

WATER = frozenset({"HOH", "DOD", "WAT", "H2O"})
COFACTORS = frozenset({
    "HEM", "HEC", "HEA", "HEB", "FAD", "FMN", "NAD", "NAP", "NDP", "NAI", "ATP", "ADP", "AMP",
    "ANP", "GTP", "GDP", "GNP", "SAM", "SAH", "COA", "ACO", "PLP", "TPP", "BCL", "CLA", "CMO",
    "OXY", "O2", "NO", "SF4", "FES", "F3S", "CU1",
})
ADDITIVES = frozenset({
    "SO4", "PO4", "GOL", "EDO", "PEG", "PG4", "PGE", "1PE", "P6G", "ACT", "FMT", "ACY", "DMS",
    "MPD", "TRS", "EPE", "MES", "CIT", "IMD", "BME", "NO3", "SCN", "NH4", "IPA", "EOH", "MSE",
    "CL", "NA", "K", "MG", "CA", "ZN", "MN", "FE", "FE2", "CU", "NI", "CO", "CD", "HG", "IOD",
    "BR", "LI", "CS", "SR", "BA", "UNX", "UNL",
})
POLAR = frozenset({"N", "O"})
APOLAR = frozenset({"C", "S"})
KIND_ORDER = ("covalent", "hbond", "hydrophobic", "contact")


@dataclass
class Ligand:
    code: str                 # código de 3 letras del PDB ("HEM", "5L4"...)
    chain: str
    resi: int
    atoms: list[Atom]
    kind: str                 # 'drug', 'cofactor' o 'additive'

    @property
    def heavy_atoms(self) -> list[Atom]:
        return [a for a in self.atoms if a.element != "H"]

    @property
    def key(self) -> str:
        return f"{self.code}/{self.chain}{self.resi}"

    def center(self) -> tuple[float, float, float]:
        atoms = self.heavy_atoms or self.atoms
        n = len(atoms)
        return (sum(a.x for a in atoms) / n, sum(a.y for a in atoms) / n, sum(a.z for a in atoms) / n)


@dataclass(frozen=True)
class Contact:
    ligand_atom: Atom
    protein_atom: Atom
    distance: float           # Å
    kind: str                 # ver KIND_ORDER


@dataclass(frozen=True)
class PocketResidue:
    """Un aminoácido del bolsillo: con qué tipo de contacto toca el fármaco y a qué distancia."""

    chain: str
    resi: int
    resn: str
    kind: str                 # el contacto más fuerte (covalent > hbond > hydrophobic > contact)
    distance: float           # Å, la del contacto más corto de ESE tipo
    contacts: int             # cuántos pares de átomos se tocan


def classify(code: str, heavy_atoms: int) -> str | None:
    code = code.upper()
    if code in WATER:
        return None
    if code in COFACTORS:
        return "cofactor"
    if code in ADDITIVES or heavy_atoms < MIN_DRUG_ATOMS:
        return "additive"
    return "drug"


def find_ligands(structure: Structure, include_additives: bool = False) -> list[Ligand]:
    """Las moléculas HETATM de la estructura (sin agua), agrupadas por residuo, en el orden
    del archivo: primero los fármacos, después los cofactores y, si se piden, los aditivos."""
    groups: dict[tuple, list[Atom]] = {}
    for atom in structure.atoms:
        if atom.hetero:
            groups.setdefault((atom.resn, atom.chain, atom.resi, atom.icode), []).append(atom)
    found = []
    for (code, chain, resi, _), atoms in groups.items():
        kind = classify(code, sum(1 for a in atoms if a.element != "H"))
        if kind is None or (kind == "additive" and not include_additives):
            continue
        found.append(Ligand(code, chain, resi, atoms, kind))
    order = {"drug": 0, "cofactor": 1, "additive": 2}
    return sorted(found, key=lambda lig: order[lig.kind])


def _kind(a: Atom, b: Atom, distance: float) -> str:
    if distance < COVALENT_CUTOFF:
        return "covalent"
    if a.element in POLAR and b.element in POLAR and distance <= HBOND_CUTOFF:
        return "hbond"
    if a.element in APOLAR and b.element in APOLAR:
        return "hydrophobic"
    return "contact"


def contacts(structure: Structure, ligand: Ligand, cutoff: float = CONTACT_CUTOFF) -> list[Contact]:
    """Pares átomo del ligando - átomo de la proteína a menos de `cutoff` Å, del más cercano al
    más lejano."""
    protein = [a for a in structure.atoms if not a.hetero and a.element != "H"]
    out = []
    for la in ligand.heavy_atoms:
        for pa in protein:
            # Descarte rápido por coordenadas antes de calcular la raíz
            if abs(la.x - pa.x) > cutoff or abs(la.y - pa.y) > cutoff or abs(la.z - pa.z) > cutoff:
                continue
            d = math.dist(la.xyz, pa.xyz)
            if d <= cutoff:
                out.append(Contact(la, pa, round(d, 2), _kind(la, pa, d)))
    return sorted(out, key=lambda c: c.distance)


def pocket(structure: Structure, ligand: Ligand, cutoff: float = CONTACT_CUTOFF) -> list[PocketResidue]:
    """Los aminoácidos que rodean al ligando, del contacto más fuerte al más débil."""
    by_residue: dict[tuple, list[Contact]] = {}
    for c in contacts(structure, ligand, cutoff):
        p = c.protein_atom
        by_residue.setdefault((p.chain, p.resi, p.resn), []).append(c)
    rank = {kind: i for i, kind in enumerate(KIND_ORDER)}
    residues = []
    for (chain, resi, resn), found in by_residue.items():
        kind = min((c.kind for c in found), key=rank.__getitem__)
        distance = min(c.distance for c in found if c.kind == kind)
        residues.append(PocketResidue(chain, resi, resn, kind, distance, len(found)))
    return sorted(residues, key=lambda r: (rank[r.kind], r.distance))
