"""
Genera dnainvent/scene/templates/b_dna.json: plantillas atómicas de un par de bases
en ADN forma B (con hidrógenos) + la transformación helicoidal entre pares.

Solo hace falta ejecutarlo una vez (el JSON ya viene incluido en el proyecto).
Requiere PyMOL open-source:  pip install pymol-open-source

La idea:
  - Para cada base X construimos la secuencia "AXA" en doble hélice.
  - Tomamos el par central (X en la hebra 1 + su complementaria en la hebra 2),
    así sus enlaces e hidrógenos son los de un nucleótido "interior".
  - Deshacemos un paso helicoidal para dejarlo en el marco del par 0.
  - El programa luego apila pares aplicando el paso k veces.
"""
import json
from pathlib import Path

import numpy as np
from pymol import cmd

OUT = Path(__file__).resolve().parent.parent / "dnainvent" / "scene" / "templates" / "b_dna.json"
COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}
BACKBONE = ["P", "O1P", "O2P", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'"]


def kabsch(src, dst):
    """Rotación R y traslación t tales que dst ≈ R @ src + t."""
    cs, cd = src.mean(0), dst.mean(0)
    h = (src - cs).T @ (dst - cd)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1, 1, d]) @ u.T
    return r, cd - r @ cs


def residues(obj):
    model = cmd.get_model(obj)
    out = {}
    for a in model.atom:
        out.setdefault((a.chain, int(a.resi)), []).append(a)
    return out, model


def strip_phosphate_h(obj):
    """PyMOL protona el fosfato; en el ADN real (pH fisiológico) está cargado (-)."""
    cmd.remove(f"{obj} and hydro and (neighbor (name O1P or name O2P))")


def main():
    # 1) Paso helicoidal, medido sobre el esqueleto de nucleótidos consecutivos
    cmd.fnab("AAAAA", name="ref", mode="DNA", form="B", dbl_helix=1)
    res, _ = residues("ref")
    a2 = {a.name: a.coord for a in res[("A", 2)]}
    a3 = {a.name: a.coord for a in res[("A", 3)]}
    src = np.array([a2[n] for n in BACKBONE])
    dst = np.array([a3[n] for n in BACKBONE])
    R, t = kabsch(src, dst)
    cmd.delete("ref")

    # El eje helicoidal es el eje de la rotación R
    w, v = np.linalg.eig(R)
    axis = np.real(v[:, np.argmin(np.abs(w - 1))])
    axis /= np.linalg.norm(axis)
    if axis @ t < 0:
        axis = -axis
    twist = np.degrees(np.arccos((np.trace(R) - 1) / 2))
    rise = float(axis @ t)
    print(f"Paso helicoidal: giro {twist:.2f}°, subida {rise:.3f} Å, eje {axis.round(3)}")

    Rinv, tinv = R.T, -R.T @ t

    templates = {}
    for base in "ATGC":
        obj = f"t{base}"
        cmd.fnab(f"A{base}A", name=obj, mode="DNA", form="B", dbl_helix=1)
        cmd.h_add(obj)
        strip_phosphate_h(obj)
        res, model = residues(obj)
        idx_of = {id(a): i for i, a in enumerate(model.atom)}

        atoms = []
        index_map = {}
        for strand, key in ((1, ("A", 2)), (2, ("B", -2))):
            for a in res[key]:
                p = np.array(a.coord) @ Rinv.T + tinv   # deshacer un paso
                index_map[idx_of[id(a)]] = len(atoms)
                part = ("phosphate" if a.name in ("P", "O1P", "O2P") else
                        "sugar" if "'" in a.name else "base")
                atoms.append({
                    "el": a.symbol, "name": a.name, "strand": strand,
                    "resn": a.resn.strip()[-1], "part": part,
                    "xyz": [round(float(c), 4) for c in p],
                })
        # Partes de hidrógenos: heredan la del átomo pesado al que se unen
        bonds = []
        for b in model.bond:
            i, j = b.index
            if i in index_map and j in index_map:
                bonds.append([index_map[i], index_map[j]])
        for i, j in bonds:
            ai, aj = atoms[i], atoms[j]
            if ai["el"] == "H":
                ai["part"] = aj["part"]
            if aj["el"] == "H":
                aj["part"] = ai["part"]
        templates[base + COMPLEMENT[base]] = {"atoms": atoms, "bonds": bonds}
        n_h = sum(a["el"] == "H" for a in atoms)
        print(f"Par {base}-{COMPLEMENT[base]}: {len(atoms)} átomos ({n_h} H), {len(bonds)} enlaces")
        cmd.delete(obj)

    data = {
        "source": "PyMOL fnab, ADN forma B, hidrógenos añadidos, fosfato desprotonado",
        "units": "angstrom",
        "step": {"R": R.round(6).tolist(), "t": t.round(5).tolist(),
                 "twist_deg": round(float(twist), 3), "rise_A": round(rise, 4),
                 "axis": axis.round(6).tolist()},
        # Enlace O3'(par k, hebra 1) -> P(par k+1, hebra 1); en la hebra 2 al revés
        "pairs": templates,
    }
    OUT.write_text(json.dumps(data, separators=(",", ":")))
    print("Guardado en", OUT)


if __name__ == "__main__":
    main()
