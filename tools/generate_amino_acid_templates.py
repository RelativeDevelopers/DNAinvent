"""
Genera dnainvent/protein/templates/amino_acids.json: dónde va cada hidrógeno de los 20
aminoácidos, en coordenadas internas, y los enlaces de cada aminoácido.

Solo hace falta ejecutarlo una vez (el JSON ya viene incluido en el proyecto).
Requiere PyMOL open-source:  pip install pymol-open-source

La idea:
  - Para cada aminoácido X se construye el tripéptido G-X-G y PyMOL le añade hidrógenos
    (pH 7: Lys y Arg con carga +, Asp y Glu con carga -, His con el H en Nε).
  - Del residuo central se guarda, para cada hidrógeno H unido al átomo pesado P:
        refs = [R2, R3]   dos átomos pesados de referencia
        r = |P-H|,  angle = R2-P-H,  dihedral = R3-R2-P-H
    con R2 y R3 unidos a P si P tiene dos vecinos pesados (posición fija), o R2 unido a P
    y R3 unido a R2 si P solo tiene uno (grupo que gira: -CH3, -OH, -NH3+, -SH). En ese
    caso se conserva la posición alternada (staggered) que elige PyMOL.
  - "-C" es el C del residuo anterior (el H del enlace peptídico va en su plano).
  - Los extremos de la cadena se calculan aparte: PyMOL los deja neutros (-NH2), pero a
    pH 7 el extremo N está cargado (-NH3+).
"""
import json
import math
import sys
from pathlib import Path

from pymol import cmd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dnainvent.constants import AMINO_ACID_CODES
from dnainvent.protein.geometry import add, angle, cross, dihedral, distance, scale, sub, unit
from dnainvent.protein.pdb import normalize_name

OUT = (Path(__file__).resolve().parent.parent / "dnainvent" / "protein" / "templates"
       / "amino_acids.json")
TETRAHEDRAL = 109.47
DIGITS = 4


def central_residue(letter: str):
    cmd.delete("all")
    cmd.fab(f"G{letter}G", "pep", ss=1)
    cmd.h_add("pep")
    model = cmd.get_model("pep")
    atoms = {i: a for i, a in enumerate(model.atom)}
    name = {i: ("-" if a.resi == "1" else "+" if a.resi == "3" else "") + normalize_name(a.name)
            for i, a in atoms.items()}
    neighbors = {i: [] for i in atoms}
    for bond in model.bond:
        i, j = bond.index
        neighbors[i].append(j)
        neighbors[j].append(i)
    return atoms, name, neighbors


GREEK = "ABGDEZH"     # α, β, γ, δ, ε, ζ, η: distancia al carbono α en la cadena lateral


def chain_rank(atom_name: str, parent: str) -> int:
    """Orden a lo largo de la cadena: N antes que C en el carbono α; en la cadena lateral,
    por la letra griega (CB < CG < CD...). El N de la prolina va el último (cierra el anillo)."""
    if parent == "CA":
        return {"N": 0, "C": 1}.get(atom_name, 2)
    if atom_name in ("N", "C"):
        return 99
    return GREEK.index(atom_name[1]) if atom_name[1:2] in GREEK else 50


def name_methylene_hydrogens(atoms, name, neighbors, mine, heavy_neighbors):
    """Nombres IUPAC de los dos H de cada grupo CH2 (HB2/HB3, HG2/HG3...; en Gly HA2/HA3).

    PyMOL los numera sin seguir la convención del PDB. La regla (comprobada con la
    ubiquitina 1D3Z, una estructura de RMN con hidrógenos): mirando el CH2 desde el vecino
    más cercano al carbono α ("anterior") hacia el siguiente, el H con diedro
    siguiente-anterior-C-H positivo es el "2" y el otro el "3".
    """
    for parent in (i for i in mine if atoms[i].symbol != "H"):
        hs = [j for j in neighbors[parent] if atoms[j].symbol == "H"]
        near = heavy_neighbors(parent)
        if len(hs) != 2 or len(near) < 2:
            continue
        pname = name[parent]
        previous, following = sorted(near, key=lambda j: chain_rank(name[j], pname))[:2]
        base = "H" + pname[1:]                      # CB -> HB, CG1 -> HG1, CA -> HA
        for h in hs:
            d = dihedral(*(tuple(atoms[k].coord) for k in (following, previous, parent, h)))
            name[h] = base + ("2" if d > 0 else "3")


def build(letter: str) -> dict:
    atoms, name, neighbors = central_residue(letter)
    mine = [i for i, a in atoms.items() if a.resi == "2"]
    heavy = [i for i in mine if atoms[i].symbol != "H"]
    order = {i: k for k, i in enumerate(sorted(atoms, key=lambda i: (atoms[i].resi != "2", i)))}

    def heavy_neighbors(i, exclude=()):
        out = [j for j in neighbors[i] if atoms[j].symbol != "H" and j not in exclude
               and (atoms[j].resi == "2" or name[j] == "-C")]
        return sorted(out, key=order.get)

    name_methylene_hydrogens(atoms, name, neighbors, mine, heavy_neighbors)
    hydrogens = []
    for h in (i for i in mine if atoms[i].symbol == "H"):
        (parent,) = neighbors[h]
        near = heavy_neighbors(parent)
        rotates = len(near) < 2
        if not rotates:
            r2, r3 = near[0], near[1]
        else:
            r2 = near[0]
            r3 = heavy_neighbors(r2, exclude=(parent,))[0]
        p, H, a2, a3 = (tuple(atoms[k].coord) for k in (parent, h, r2, r3))
        hydrogens.append({
            "name": name[h], "parent": name[parent], "refs": [name[r2], name[r3]],
            "r": round(distance(p, H), DIGITS),
            "angle": round(angle(a2, p, H), DIGITS),
            "dihedral": round(dihedral(a3, a2, p, H), DIGITS),
            "rotates": rotates,
        })

    bonds = sorted({tuple(sorted((name[i], name[j]))) for i in mine for j in neighbors[i]
                    if atoms[j].resi == "2"})
    return {
        "heavy": [name[i] for i in heavy],
        "bonds": [list(b) for b in bonds],
        "hydrogens": hydrogens,
    }


def n_terminus(templates: dict) -> dict:
    """Extremo N a pH 7: -NH3+ (3 H alternados) o, en la prolina, -NH2+ (2 H)."""
    r = next(h["r"] for h in templates["ALA"]["hydrogens"] if h["name"] == "H")
    regular = [{"name": f"H{k + 1}", "parent": "N", "refs": ["CA", "C"], "r": r,
                "angle": TETRAHEDRAL, "dihedral": phi, "rotates": True} for k, phi in enumerate((60.0, 180.0, -60.0))]

    # Prolina: el N ya tiene dos vecinos (CA y CD); los dos H completan el tetraedro.
    _, name, _ = central_residue("P")
    model = cmd.get_model("pep and resi 2")
    xyz = {normalize_name(a.name): tuple(a.coord) for a in model.atom}
    n, ca, cd = xyz["N"], xyz["CA"], xyz["CD"]
    bisector = unit(add(unit(sub(ca, n)), unit(sub(cd, n))))
    normal = unit(cross(sub(ca, n), sub(cd, n)))
    half = math.radians(TETRAHEDRAL / 2)
    proline = []
    for k, sign in enumerate((1, -1)):
        direction = add(scale(bisector, -math.cos(half)), scale(normal, sign * math.sin(half)))
        h = add(n, scale(direction, r))
        proline.append({"name": f"H{k + 1}", "parent": "N", "refs": ["CA", "CD"], "r": r,
                        "angle": round(angle(ca, n, h), DIGITS),
                        "dihedral": round(dihedral(cd, ca, n, h), DIGITS), "rotates": False})
    return {"regular": regular, "proline": proline}


def main():
    templates = {}
    for letter in sorted(AMINO_ACID_CODES):
        resn = AMINO_ACID_CODES[letter].upper()
        templates[resn] = build(letter)
        t = templates[resn]
        print(f"{resn}: {len(t['heavy'])} átomos pesados, {len(t['hydrogens'])} H, "
              f"{len(t['bonds'])} enlaces")
    data = {
        "source": "PyMOL fab G-X-G + h_add (pH 7: Lys+, Arg+, Asp-, Glu-, His con H en NE2)",
        "units": {"r": "angstrom", "angle": "degrees", "dihedral": "degrees"},
        "convention": "H = place(R3, R2, parent, r, angle, dihedral); '-C' = C del residuo anterior",
        "residues": templates,
        "n_terminus": n_terminus(templates),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(data, indent=1, ensure_ascii=False), encoding="utf-8")
    print("Guardado en", OUT)


if __name__ == "__main__":
    main()
