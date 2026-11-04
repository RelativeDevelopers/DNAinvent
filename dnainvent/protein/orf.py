"""Marcos de lectura abiertos (ORF, del inglés Open Reading Frame).

Una secuencia de ADN se puede leer de 3 en 3 de seis maneras distintas (seis "pautas"):

    hebra +   empezando en la base 1, 2 o 3       ATG CGT TAA ...   /   TGC GTT AA ...
    hebra -   lo mismo en la complementaria inversa (la otra hebra, leída 5'->3')

Un ORF es un tramo que empieza en un codón de inicio (ATG) y termina en un codón de parada
(TAA, TAG o TGA) dentro de la misma pauta. El gen de verdad suele ser el ORF más largo:
en una secuencia al azar aparece un codón de parada cada ~21 codones (3 de 64), así que
los ORF largos casi nunca son casualidad.

Si un ORF llega al final de la secuencia sin encontrar parada, también se devuelve, marcado
como incompleto: puede ser un gen cortado (el usuario ha pegado solo una parte).

Las posiciones (start, end) son siempre las de la secuencia original, hebra +, desde 0 y con
end excluido, también para los ORF de la hebra -: así se pueden señalar en la doble hélice.

Aviso biológico: en organismos como el ser humano los genes están troceados (exones e
intrones). Esto funciona con la secuencia codificante (CDS o ADNc), no con el ADN genómico.
"""

from dataclasses import dataclass

from ..constants import CODON_TABLE, START_CODON, STOP_CODONS
from ..errors import InvalidParameterError, InvalidSequenceError
from ..sequences import reverse_complement
from ..validation import validate_dna

# El código genético está escrito para ARN (U); aquí se lee ADN (T).
_DNA_CODONS = {codon.replace("U", "T"): amino for codon, amino in CODON_TABLE.items()}
_DNA_START = START_CODON.replace("U", "T")
_DNA_STOPS = frozenset(codon.replace("U", "T") for codon in STOP_CODONS)


@dataclass(frozen=True)
class ORF:
    """Un marco de lectura abierto.

    strand      +1 (la secuencia tal cual) o -1 (complementaria inversa)
    frame       0, 1 o 2: en qué base empieza la lectura dentro de esa hebra
    start, end  tramo en la secuencia original (hebra +, desde 0, end excluido),
                incluido el codón de parada si lo hay
    protein     aminoácidos en código de una letra, sin la parada
    stop_codon  "TAA", "TAG" o "TGA"; None si llega al final sin parada (incompleto)
    """

    strand: int
    frame: int
    start: int
    end: int
    protein: str
    stop_codon: str | None

    @property
    def length(self) -> int:
        """Número de aminoácidos."""
        return len(self.protein)

    @property
    def is_complete(self) -> bool:
        return self.stop_codon is not None

    @property
    def label(self) -> str:
        """Nombre habitual de la pauta: "+1", "+2", "+3", "-1", "-2" o "-3"."""
        return f"{'+' if self.strand > 0 else '-'}{self.frame + 1}"


def find_orfs(dna: str, min_length: int = 1) -> list[ORF]:
    """Todos los ORF de las 6 pautas con al menos min_length aminoácidos,
    del más largo al más corto (ver _rank para los empates).

    Dentro de una pauta, un ATG que cae dentro de un ORF ya abierto no empieza otro: se
    toma el primero, que da la proteína más larga (como hace el ribosoma).
    """
    dna = validate_dna(dna)
    if min_length < 1:
        raise InvalidParameterError("min_length_too_small", minimum=1)

    orfs = []
    for strand, seq in ((1, dna), (-1, reverse_complement(dna))):
        for frame in range(3):
            orfs.extend(_scan_frame(seq, strand, frame, len(dna)))
    return sorted((orf for orf in orfs if orf.length >= min_length), key=_rank)


def best_orf(dna: str) -> ORF:
    """El ORF más probable de ser el gen: el más largo de las 6 pautas."""
    orfs = find_orfs(dna)
    if not orfs:
        raise InvalidSequenceError("no_orf")
    return orfs[0]


def _rank(orf: ORF) -> tuple:
    """Orden: más largo primero; a igual longitud, completo antes que incompleto, hebra +
    antes que -, y después por pauta y posición (para que el resultado sea siempre el mismo)."""
    return (-orf.length, not orf.is_complete, -orf.strand, orf.frame, orf.start)


def _scan_frame(seq: str, strand: int, frame: int, total: int) -> list[ORF]:
    found = []
    start = None                # posición del ATG del ORF abierto (en seq), o None
    amino_acids = []
    for position in range(frame, len(seq) - 2, 3):
        codon = seq[position:position + 3]
        if start is None:
            if codon == _DNA_START:
                start, amino_acids = position, [_DNA_CODONS[codon]]
            continue
        if codon in _DNA_STOPS:
            found.append(_make(strand, frame, start, position + 3, amino_acids, codon, total))
            start = None
        else:
            amino_acids.append(_DNA_CODONS[codon])
    if start is not None:       # llegó al final sin parada
        last = start + 3 * len(amino_acids)
        found.append(_make(strand, frame, start, last, amino_acids, None, total))
    return found


def _make(strand, frame, start, end, amino_acids, stop, total) -> ORF:
    if strand < 0:              # de la complementaria inversa a la hebra +
        start, end = total - end, total - start
    return ORF(strand, frame, start, end, "".join(amino_acids), stop)
