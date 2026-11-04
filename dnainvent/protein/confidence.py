"""Confianza de una predicción (pLDDT, de 0 a 100), en las cuatro franjas de AlphaFold.

    >= 90   very_high   muy alta: posición casi exacta, incluso las cadenas laterales
    >= 70   confident   buena: el esqueleto está bien
    >= 50   low         baja: tómalo con cautela
    <  50   very_low    muy baja: suele ser una zona desordenada (no tiene forma fija)

Los colores son los que usa AlphaFold DB, para que coincidan con su web.
"""

BANDS = (                              # (mínimo, clave, color)
    (90.0, "very_high", "#0053d6"),
    (70.0, "confident", "#65cbf3"),
    (50.0, "low", "#ffdb13"),
    (0.0, "very_low", "#ff7d45"),
)


def band(plddt: float) -> str:
    return next(key for minimum, key, _ in BANDS if plddt >= minimum)


def color(plddt: float) -> str:
    return next(c for minimum, _, c in BANDS if plddt >= minimum)


def i18n_key(plddt: float) -> str:
    return f"protein.confidence.{band(plddt)}"


def per_residue(structure) -> list[float]:
    """pLDDT de cada aminoácido (el de su carbono α)."""
    return [r.atoms["CA"].bfactor for r in structure.residues() if "CA" in r.atoms]


def mean(structure) -> float:
    values = per_residue(structure)
    return sum(values) / len(values) if values else 0.0


def counts(structure) -> dict[str, int]:
    """Aminoácidos en cada franja: {"very_high": 143, "confident": 3, ...}"""
    out = {key: 0 for _, key, _ in BANDS}
    for value in per_residue(structure):
        out[band(value)] += 1
    return out
