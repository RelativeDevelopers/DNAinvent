"""Resumen de la molécula para el visor web: todo lo que necesita para dibujarla, en tablas.

El visor de escritorio genera los puntos con NumPy. El web los genera en el navegador, pero
sin física: recibe las posiciones de los átomos y, para cada elemento (H, C, N, O, P), el
resumen de AtomicWorld (atomicworld.scene.describe): la tabla de radios de cada subnivel,
las casillas de Hund, las posiciones de los nucleones... y solo tiene que interpolar.
Así la física está en un único sitio.

    describe_structure("ATGCGTACGTTAGC")         # -> dict listo para JSON
    describe_structure(dna, start=120, count=60)  # un tramo de una secuencia larga
    describe_protein(build_protein(estructura))   # lo mismo para una proteína (scene/protein.py)
"""
from functools import lru_cache

from atomicworld.scene.describe import describe_cached

from .elements import atom_type
from .limits import RECOMMENDED_PAIRS, RECOMMENDED_RESIDUES
from .molecule import build_molecule

PART_CODES = {"base": "b", "sugar": "s", "phosphate": "p"}
PROTEIN_PART_CODES = {"backbone": "k", "sidechain": "c"}
ELEMENT_ORDER = "HCNOPS"


@lru_cache(maxsize=None)
def describe_element(symbol: str) -> dict:
    """Resumen de AtomicWorld del isótopo más abundante del elemento (~20 kB)."""
    t = atom_type(symbol)
    return describe_cached(t.number, t.isotope.a)


def describe_structure(dna: str, start: int = 0, count: int | None = None) -> dict:
    mol = build_molecule(dna, start, count)
    symbols = sorted(set(mol.elements), key=ELEMENT_ORDER.index)
    return {
        "sequence": mol.sequence,
        "start": mol.start,
        "length": mol.length,
        "recommended_pairs": RECOMMENDED_PAIRS,
        "atoms": {
            "count": mol.n_atoms,
            # 3 decimales en nm = 0,01 Å: de sobra para dibujar
            "positions_nm": [round(float(v), 3) for v in mol.positions.ravel()],
            "elements": "".join(mol.elements),        # todos los símbolos son de una letra
            "names": mol.names,
            "strand": "".join(str(s) for s in mol.strand),
            "pair": [int(k) for k in mol.pair_index],
            "resn": "".join(mol.resn),
            "part": "".join(PART_CODES[p] for p in mol.part),
        },
        "bonds": [int(v) for v in mol.bonds.ravel()],
        "hbonds": [int(v) for v in mol.hbonds.ravel()],
        "elements": {s: describe_element(s) for s in symbols},
    }


def describe_protein(mol) -> dict:
    """Una ProteinMolecule (scene/protein.py) para el visor web, con la misma forma que
    describe_structure. Lo propio de las proteínas:

        residue      nº de aminoácido de cada átomo en la cadena completa (desde 0)
        part         "k" esqueleto (backbone) o "c" cadena lateral
        resn         aminoácido en 3 letras de cada átomo
        confidence   pLDDT de cada átomo (null si no es una predicción)
    """
    symbols = sorted(set(mol.elements), key=ELEMENT_ORDER.index)
    return {
        "kind": "protein",
        "sequence": mol.sequence,
        "start": mol.start,
        "length": mol.length,
        "recommended": RECOMMENDED_RESIDUES,
        "atoms": {
            "count": mol.n_atoms,
            "positions_nm": [round(float(v), 3) for v in mol.positions.ravel()],
            "elements": "".join(mol.elements),        # H, C, N, O, S: todos de una letra
            "names": mol.names,
            "residue": [int(k) for k in mol.pair_index],
            "resn": mol.resn,
            "part": "".join(PROTEIN_PART_CODES[p] for p in mol.part),
            "confidence": (None if mol.confidence is None
                           else [round(float(c), 1) for c in mol.confidence]),
        },
        "bonds": [int(v) for v in mol.bonds.ravel()],
        "hbonds": [],
        "elements": {s: describe_element(s) for s in symbols},
    }
