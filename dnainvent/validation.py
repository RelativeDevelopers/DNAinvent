"""Limpieza y validación de secuencias. Todas las funciones de lógica pasan por aquí."""

from .constants import BITS_PER_BASE, DNA_BASES, PROTEIN_LETTERS, RNA_BASES
from .errors import InvalidSequenceError


def clean(sequence: str) -> str:
    """Quita espacios y saltos de línea y pasa a mayúsculas."""
    return "".join(sequence.split()).upper()


def _validate(sequence: str, alphabet: frozenset, kind: str) -> str:
    sequence = clean(sequence)
    if not sequence:
        raise InvalidSequenceError("empty_sequence", kind=kind)
    invalid = sorted(set(sequence) - alphabet)
    if invalid:
        raise InvalidSequenceError("invalid_characters", kind=kind,
                                   chars=invalid, allowed=sorted(alphabet))
    return sequence


def validate_dna(sequence: str) -> str:
    return _validate(sequence, DNA_BASES, "dna")


def validate_rna(sequence: str) -> str:
    return _validate(sequence, RNA_BASES, "rna")


def validate_protein(sequence: str) -> str:
    """Proteína en código de una letra (los 20 aminoácidos y los comodines B, Z y X)."""
    return _validate(sequence, PROTEIN_LETTERS, "protein")


def validate_bits(bits: str) -> str:
    """Cadena de bits que se pueda leer de 2 en 2 como bases."""
    bits = _validate(bits, frozenset("01"), "binary")
    if len(bits) % BITS_PER_BASE:
        raise InvalidSequenceError("odd_bit_length", length=len(bits), multiple=BITS_PER_BASE)
    return bits
