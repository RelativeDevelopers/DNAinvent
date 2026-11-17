"""La mutación en la proteína 3D (fase 1e): dónde cae el aminoácido que cambia, en la estructura
de AlphaFold de la proteína sana.

    site = protein3d.variant_site(conn, variant_id, AlphaFoldDB())
    site["label"]       # 'Glu7Val'
    site["burial"]      # 'exposed': en la superficie de la proteína
    site["changes"]     # ['loses_negative', 'polar_to_hydrophobic']

Devuelve datos, sin textos traducidos: la consola los pinta y la web (fase 2) hará lo mismo.

Lo que NO hace, a propósito: predecir la forma de la proteína con la mutación. AlphaFold y
ESMFold no captan bien el efecto de un solo aminoácido (casi siempre dan la misma forma) y eso
haría creer que la mutación no hace nada. Se enseña la proteína sana con el sitio marcado y se
explica qué cambia en el aminoácido y dónde está.

Numeración: la de la proteína de RefSeq de la variante (HGVS, desde la metionina inicial). La
estructura es la de UniProt; si las dos secuencias no son iguales se alinean, y si el aminoácido
de la estructura no es el que dice la variante, no se marca nada (mejor nada que el sitio
equivocado).

AlphaFold DB no publica modelos de proteínas de más de 2.700 aminoácidos (solo fragmentos en
sus descargas): BRCA2, HTT y DMD se quedan sin vista 3D.
"""

import re
import sqlite3

from ..constants import AMINO_ACID_CODES
from ..errors import DNAinventError
from ..protein import confidence
from ..protein.compare import pairs_by_alignment
from ..protein.geometry import distance
from ..protein.pdb import Structure
from ..scene.limits import RECOMMENDED_RESIDUES
from .variants import locate

MAX_ALPHAFOLD_LENGTH = 2700          # aminoácidos: por encima, AlphaFold DB no tiene el modelo
NEIGHBOR_RADIUS = 10.0               # Å entre carbonos α, para saber si está dentro o fuera
BURIED, EXPOSED = 20, 12             # vecinos: ≥ 20 dentro, ≤ 12 en la superficie
CONTACT_RADIUS = 4.5                 # Å entre átomos, para los aminoácidos que lo tocan
CONTACTS_SHOWN = 6

THREE_TO_ONE = {v: k for k, v in AMINO_ACID_CODES.items()}
HGVS_P = re.compile(r"p\.\(?([A-Z][a-z]{2})(\d+)(.*?)\)?$")

# Propiedades de los aminoácidos (las que se explican en la ficha)
CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1}
HYDROPHOBIC = set("AVLIMFWC")
VOLUME = {"G": 60, "A": 89, "S": 89, "C": 109, "D": 111, "P": 113, "N": 114, "T": 116,
          "E": 138, "V": 140, "Q": 144, "H": 153, "M": 163, "I": 167, "L": 167, "K": 169,
          "R": 174, "F": 190, "Y": 194, "W": 228}          # Å³ (Zamyatnin, 1972)
SIZE_CHANGE = 40                                         # Å³ para decir «más grande / pequeño»


# ---------------------------------------------------------------------- el cambio
def variant_change(conn: sqlite3.Connection, variant_id: int) -> dict | None:
    """Qué aminoácido cambia: {'position', 'ref', 'alt', 'kind', 'label'}, o None si la variante
    no tiene un aminoácido concreto (p.?, fuera del CDS...).

    Primero con el codón calculado (sustituciones de una letra); si no, con el HGVS de la
    proteína que da ClinVar (deleciones, desplazamientos del marco...)."""
    try:
        effect, _, _ = locate(conn, variant_id)
    except DNAinventError:
        effect = None
    if effect is not None:
        if effect.kind == "synonymous":
            return None
        label = effect.hgvs_p[2:]
        return {"position": effect.codon_number, "ref": effect.ref_aa,
                "alt": effect.alt_aa, "kind": effect.kind, "label": label}
    row = conn.execute("SELECT hgvs_p FROM variants WHERE id = ?", (variant_id,)).fetchone()
    return parse_hgvs_p(row["hgvs_p"] if row else None)


def parse_hgvs_p(hgvs_p: str | None) -> dict | None:
    """'p.Phe508del' -> {'position': 508, 'ref': 'F', 'alt': None, 'kind': 'deletion', ...}"""
    match = HGVS_P.match(hgvs_p or "")
    if not match or match.group(1) not in THREE_TO_ONE and match.group(1) != "Ter":
        return None
    ref = "*" if match.group(1) == "Ter" else THREE_TO_ONE[match.group(1)]
    rest = match.group(3)
    alt, kind = None, "other"
    if "fs" in rest:
        kind = "frameshift"
    elif rest.startswith("del") or rest.startswith("_") and "del" in rest:
        kind = "deletion"
    elif rest.startswith("dup") or "ins" in rest:
        kind = "insertion"
    elif rest in ("Ter", "*"):
        kind, alt = "nonsense", "*"
    elif rest in THREE_TO_ONE:
        alt = THREE_TO_ONE[rest]
        kind = "stop_lost" if ref == "*" else "missense"
    elif rest == "=":
        return None
    return {"position": int(match.group(2)), "ref": ref, "alt": alt, "kind": kind,
            "label": hgvs_p.replace("(", "").replace(")", "")[2:]}


def change_properties(ref: str, alt: str | None) -> list[str]:
    """Claves de lo que cambia del aminoácido (para i18n: cards.protein3d.changes.*)."""
    if not alt or alt == "*" or ref == "*" or ref not in VOLUME or alt not in VOLUME:
        return []
    out = []
    before, after = CHARGE.get(ref, 0), CHARGE.get(alt, 0)
    if before != after:
        if before and after:
            out.append("charge_reversed")
        elif before < 0:
            out.append("loses_negative")
        elif before > 0:
            out.append("loses_positive")
        elif after < 0:
            out.append("gains_negative")
        else:
            out.append("gains_positive")
    if ref in HYDROPHOBIC and alt not in HYDROPHOBIC:
        out.append("hydrophobic_to_polar")
    elif ref not in HYDROPHOBIC and alt in HYDROPHOBIC:
        out.append("polar_to_hydrophobic")
    if VOLUME[alt] - VOLUME[ref] >= SIZE_CHANGE:
        out.append("bigger")
    elif VOLUME[ref] - VOLUME[alt] >= SIZE_CHANGE:
        out.append("smaller")
    if alt == "P":
        out.append("to_proline")
    if ref == "G":
        out.append("from_glycine")
    if "C" in (ref, alt):
        out.append("cysteine")
    return out


# ---------------------------------------------------------------------- en la estructura
def map_position(refseq: str, structure_sequence: str, position: int) -> int | None:
    """Posición de la variante (desde 1, RefSeq) -> índice del aminoácido en la estructura
    (desde 0), o None si no tiene pareja."""
    if refseq == structure_sequence:
        return position - 1 if 0 < position <= len(refseq) else None
    pairs = dict(pairs_by_alignment(refseq, structure_sequence))
    return pairs.get(position - 1)


def burial(structure: Structure, index: int) -> tuple[str, int]:
    """('buried' | 'partial' | 'exposed', vecinos): cuántos carbonos α hay a menos de 10 Å."""
    residues = structure.residues()
    ca = residues[index].atoms.get("CA")
    if ca is None:
        return "partial", 0
    neighbors = sum(1 for k, r in enumerate(residues) if k != index and "CA" in r.atoms
                    and distance(r.atoms["CA"].xyz, ca.xyz) < NEIGHBOR_RADIUS)
    kind = "buried" if neighbors >= BURIED else "exposed" if neighbors <= EXPOSED else "partial"
    return kind, neighbors


def contacts(structure: Structure, index: int) -> list[dict]:
    """Aminoácidos que tocan al de la variante (algún átomo a menos de 4,5 Å), sin contar los
    dos de al lado en la cadena. Los más cercanos primero."""
    residues = structure.residues()
    own = list(residues[index].atoms.values())
    found = []
    for k, residue in enumerate(residues):
        if abs(k - index) <= 1:
            continue
        closest = min(distance(a.xyz, b.xyz) for a in own for b in residue.atoms.values())
        if closest < CONTACT_RADIUS:
            found.append({"index": k, "number": residue.resi, "residue": residue.letter,
                          "distance": round(closest, 1)})
    return sorted(found, key=lambda c: c["distance"])


def analyze(structure: Structure, refseq: str, change: dict) -> dict:
    """Todo lo que la ficha dice del sitio. 'available' False con 'reason' si no se puede."""
    residues = structure.residues()
    index = map_position(refseq, structure.sequence(), change["position"])
    if index is None or index >= len(residues):
        return {"available": False, "reason": "not_in_structure"}
    residue = residues[index]
    if change["ref"] != "*" and residue.letter != change["ref"]:
        return {"available": False, "reason": "mismatch"}
    plddt = residue.atoms["CA"].bfactor if "CA" in residue.atoms and structure.is_prediction else None
    kind, neighbors = burial(structure, index)
    if plddt is not None and plddt < 50:
        kind = "disordered"                    # sin forma fija: dentro/fuera no tiene sentido
    length = len(residues)
    count = min(RECOMMENDED_RESIDUES, length)
    start = min(max(0, index - count // 2), length - count)
    site = {
        "available": True, "index": index, "number": residue.resi, "length": length,
        "label": change["label"], "kind": change["kind"], "ref": change["ref"], "alt": change["alt"],
        "plddt": plddt, "confidence": confidence.band(plddt) if plddt is not None else None,
        "burial": kind, "neighbors": neighbors,
        "contacts": contacts(structure, index)[:CONTACTS_SHOWN],
        "changes": change_properties(change["ref"], change["alt"]),
        "start": start, "count": count,
    }
    if change["kind"] in ("nonsense", "frameshift"):
        site["lost"] = length - index         # aminoácidos que ya no se fabrican (o cambian)
    return site


def variant_site(conn: sqlite3.Connection, variant_id: int, afdb) -> dict:
    """El sitio de una variante en la estructura de AlphaFold de su proteína.

    Siempre devuelve un diccionario: con 'available' False y 'reason' ('no_change',
    'no_uniprot', 'too_long', 'service', 'not_in_structure', 'mismatch') si no se puede."""
    change = variant_change(conn, variant_id)
    if change is None:
        return {"available": False, "reason": "no_change"}
    gene = conn.execute(
        """SELECT g.symbol, g.uniprot_id, s.residues FROM variants v JOIN genes g ON g.id = v.gene_id
           LEFT JOIN sequences s ON s.gene_id = g.id AND s.mol_type = 'protein'
           WHERE v.id = ?""", (variant_id,)).fetchone()
    base = {"gene": gene["symbol"], "uniprot": gene["uniprot_id"], "label": change["label"]}
    if not gene["uniprot_id"] or not gene["residues"]:
        return {**base, "available": False, "reason": "no_uniprot"}
    if len(gene["residues"]) > MAX_ALPHAFOLD_LENGTH:
        return {**base, "available": False, "reason": "too_long", "length": len(gene["residues"])}
    try:
        entry, structure = afdb.structure(gene["uniprot_id"])
    except DNAinventError as exc:
        return {**base, "available": False, "reason": "service", "error": exc}
    site = analyze(structure, gene["residues"], change)
    return {**base, **site, "model_id": entry.model_id, "version": entry.version,
            "structure": structure}
