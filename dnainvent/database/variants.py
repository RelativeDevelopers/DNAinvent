"""Situar una variante en su ARNm y compararla con la secuencia sana.

Ejemplo, la anemia falciforme (HbS): NM_000518.5:c.20A>T
  * "c.20" es la base 20 contando desde la A del ATG inicial del CDS.
  * En el ARNm, eso es la posición cds_start + 20 - 1 (70 en NM_000518.5, cuyo CDS empieza en 51).
  * La secuencia de referencia en esa posición TIENE que ser A. Si no lo es, la variante y la
    secuencia no corresponden (otra versión del ARNm) y se avisa en vez de inventar.
  * El codón 7 pasa de GAG a GTG: ácido glutámico -> valina (p.Glu7Val).

Por ahora solo sustituciones de una base dentro del CDS (c.N X>Y), que son la mayoría de las
variantes de ClinVar. Usa el código genético estándar (genes nucleares).

Numeración: HGVS cuenta la metionina inicial como aminoácido 1. En textos antiguos HbS aparece
como Glu6Val porque no la contaban (en la hemoglobina madura se elimina): es la misma mutación.
"""

import re
import sqlite3
from dataclasses import dataclass

from ..constants import AMINO_ACID_CODES, CODON_TABLE
from ..errors import InvalidParameterError, ServiceError

SIMPLE_SNV = re.compile(r"c\.(\d+)([ACGT])>([ACGT])")
CONTEXT_CODONS = 5          # codones a cada lado en la comparación


def amino_acid(codon: str) -> str:
    """Codón de ADN -> letra del aminoácido ('*' = parada)."""
    letter = CODON_TABLE[codon.replace("T", "U")]
    return "*" if letter == "_" else letter


def three_letter(letter: str) -> str:
    return "Ter" if letter == "*" else AMINO_ACID_CODES[letter]


@dataclass(frozen=True)
class VariantEffect:
    accession: str            # "NM_000518.5"
    hgvs_c: str               # "c.20A>T"
    position: int             # en el ARNm, desde 1 (70)
    codon_number: int         # 7
    ref_codon: str            # "GAG"
    alt_codon: str            # "GTG"
    ref_aa: str               # "E"
    alt_aa: str               # "V"
    kind: str                 # 'synonymous', 'missense', 'nonsense' o 'stop_lost'

    @property
    def hgvs_p(self) -> str:
        alt = "=" if self.kind == "synonymous" else three_letter(self.alt_aa)
        return f"p.{three_letter(self.ref_aa)}{self.codon_number}{alt}"


@dataclass(frozen=True)
class SideBySide:
    """Un tramo del CDS, sano y con la variante, codón a codón."""

    first_codon: int          # número del primer codón del tramo
    healthy: list[str]        # codones
    variant: list[str]
    healthy_protein: str      # una letra por codón
    variant_protein: str
    changed_codon: int        # índice (en las listas) del codón que cambia
    changed_base: int         # posición de la base cambiada dentro de ese codón (0-2)


def locate(conn: sqlite3.Connection, variant_id: int) -> tuple[VariantEffect, str, int]:
    """Calcula el efecto de una sustitución c.N X>Y comprobándola contra la secuencia.

    Devuelve (efecto, secuencia del ARNm, cds_start). Guarda la posición en variants.position.
    La secuencia tiene que estar descargada (NCBIProvider.ensure_residues)."""
    row = conn.execute(
        """SELECT v.hgvs_c, s.residues, s.cds_start, s.cds_end, s.accession, s.version
           FROM variants v JOIN sequences s ON s.id = v.sequence_id WHERE v.id = ?""",
        (variant_id,)).fetchone()
    if row is None:
        raise InvalidParameterError("variant_not_found", variant=variant_id)
    hgvs = row["hgvs_c"] or ""
    accession = f"{row['accession']}.{row['version']}"
    match = SIMPLE_SNV.fullmatch(hgvs)
    if not match:
        raise InvalidParameterError("variant_not_supported", hgvs=hgvs)
    if not row["residues"] or row["cds_start"] is None:
        raise ServiceError("sequence_not_found", accession=accession)

    c_pos, ref, alt = int(match.group(1)), match.group(2), match.group(3)
    residues, cds_start, cds_end = row["residues"], row["cds_start"], row["cds_end"]
    position = cds_start + c_pos - 1
    if position > cds_end:
        raise InvalidParameterError("variant_outside_cds", hgvs=hgvs, start=cds_start, end=cds_end)
    found = residues[position - 1]
    if found != ref:
        raise InvalidParameterError("variant_mismatch", hgvs=hgvs, accession=accession,
                                    position=position, found=found, expected=ref)

    codon_number = (c_pos - 1) // 3 + 1
    codon_start = cds_start + (codon_number - 1) * 3
    ref_codon = residues[codon_start - 1:codon_start + 2]
    offset = position - codon_start
    alt_codon = ref_codon[:offset] + alt + ref_codon[offset + 1:]
    ref_aa, alt_aa = amino_acid(ref_codon), amino_acid(alt_codon)
    kind = ("synonymous" if ref_aa == alt_aa else "nonsense" if alt_aa == "*"
            else "stop_lost" if ref_aa == "*" else "missense")

    with conn:
        conn.execute("UPDATE variants SET position = ?, ref_allele = ?, alt_allele = ? WHERE id = ?",
                     (position, ref, alt, variant_id))
    effect = VariantEffect(accession, hgvs, position, codon_number, ref_codon, alt_codon,
                           ref_aa, alt_aa, kind)
    return effect, residues, cds_start


def side_by_side(effect: VariantEffect, residues: str, cds_start: int,
                 context: int = CONTEXT_CODONS) -> SideBySide:
    """Los codones alrededor del cambio, sanos y con la variante."""
    first = max(1, effect.codon_number - context)
    healthy = []
    for number in range(first, effect.codon_number + context + 1):
        start = cds_start - 1 + (number - 1) * 3
        codon = residues[start:start + 3]
        if len(codon) < 3:
            break
        healthy.append(codon)
    index = effect.codon_number - first
    variant = list(healthy)
    variant[index] = effect.alt_codon
    return SideBySide(
        first_codon=first, healthy=healthy, variant=variant,
        healthy_protein="".join(amino_acid(c) for c in healthy),
        variant_protein="".join(amino_acid(c) for c in variant),
        changed_codon=index,
        changed_base=effect.position - (cds_start + (effect.codon_number - 1) * 3),
    )
