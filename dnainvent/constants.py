"""Tablas de datos biológicos. Solo datos, sin lógica."""

DNA_BASES = frozenset("ACGT")
RNA_BASES = frozenset("ACGU")

# Pares de bases de Watson-Crick.
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

START_CODON = "AUG"
STOP_CODONS = frozenset({"UAA", "UAG", "UGA"})

# Código genético estándar: codón de ARN -> aminoácido (una letra). '_' = parada.
CODON_TABLE = {
    "AUA": "I", "AUC": "I", "AUU": "I", "AUG": "M",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACU": "T",
    "AAC": "N", "AAU": "N", "AAA": "K", "AAG": "K",
    "AGC": "S", "AGU": "S", "AGA": "R", "AGG": "R",
    "CUA": "L", "CUC": "L", "CUG": "L", "CUU": "L",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCU": "P",
    "CAC": "H", "CAU": "H", "CAA": "Q", "CAG": "Q",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGU": "R",
    "GUA": "V", "GUC": "V", "GUG": "V", "GUU": "V",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCU": "A",
    "GAC": "D", "GAU": "D", "GAA": "E", "GAG": "E",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGU": "G",
    "UCA": "S", "UCC": "S", "UCG": "S", "UCU": "S",
    "UUC": "F", "UUU": "F", "UUA": "L", "UUG": "L",
    "UAC": "Y", "UAU": "Y", "UAA": "_", "UAG": "_",
    "UGC": "C", "UGU": "C", "UGA": "_", "UGG": "W",
}

# Códigos de 3 letras (internacionales). Los nombres completos están en i18n/locales/.
AMINO_ACID_CODES = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
}

# Codificación de 2 bits por base (4 bases = 2^2 combinaciones).
# El orden no es casual: cada base y su complementaria son opuestas bit a bit
# (A=00 <-> T=11, C=01 <-> G=10), así que invertir los bits da la complementaria.
BASE_TO_BINARY = {"A": "00", "C": "01", "G": "10", "T": "11"}
BINARY_TO_BASE = {bits: base for base, bits in BASE_TO_BINARY.items()}
BITS_PER_BASE = 2


# ---------- Proteínas ----------

# Tipos de aminoácido según su cadena lateral (lo que decide dónde se colocan al plegarse):
#   hydrophobic  huyen del agua: suelen quedar en el interior de la proteína
#   special      glicina (la más flexible) y prolina (rígida, rompe hélices)
#   polar        sin carga, pero atraen el agua
#   positive     carga + a pH 7
#   negative     carga - a pH 7
AMINO_ACID_CLASSES = {
    **dict.fromkeys("AVLIMFW", "hydrophobic"),
    **dict.fromkeys("GP", "special"),
    **dict.fromkeys("STNQYC", "polar"),
    **dict.fromkeys("KRH", "positive"),
    **dict.fromkeys("DE", "negative"),
}

# Letras que aceptan las funciones de proteínas: los 20 aminoácidos estándar y los tres
# comodines que usa UniProt: B (D o N), Z (E o Q) y X (desconocido).
PROTEIN_LETTERS = frozenset(AMINO_ACID_CODES) | {"B", "Z", "X"}

# BLOSUM62 (Henikoff y Henikoff, 1992): cuánto suma o resta encontrar el aminoácido de la
# fila frente al de la columna en un alineamiento. Positivo = cambio frecuente en la
# evolución (I <-> V: +3); negativo = raro (W <-> D: -4). Es la matriz que usa BLAST.
# Valores comprobados contra Bio.Align.substitution_matrices.load("BLOSUM62").
BLOSUM62_ORDER = "ARNDCQEGHILKMFPSTWYVBZX"
_BLOSUM62_ROWS = {
    "A": " 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0",
    "R": "-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1",
    "N": "-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1",
    "D": "-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1",
    "C": " 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2",
    "Q": "-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1",
    "E": "-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1",
    "G": " 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1",
    "H": "-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1",
    "I": "-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1",
    "L": "-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1",
    "K": "-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1",
    "M": "-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1",
    "F": "-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1",
    "P": "-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2",
    "S": " 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0",
    "T": " 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0",
    "W": "-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2",
    "Y": "-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1",
    "V": " 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1",
    "B": "-2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1",
    "Z": "-1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1",
    "X": " 0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1",
}
# BLOSUM62["W"]["W"] == 11
BLOSUM62 = {row: dict(zip(BLOSUM62_ORDER, map(int, values.split())))
            for row, values in _BLOSUM62_ROWS.items()}
