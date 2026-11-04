"""Tests de la fase 1 de proteínas: alineamiento de proteínas, ORF y clasificación (sin red)."""

import json
import random
from functools import lru_cache

import pytest

from dnainvent import InvalidParameterError, InvalidSequenceError, alignment, sequences
from dnainvent.constants import BLOSUM62, BLOSUM62_ORDER, CODON_TABLE
from dnainvent.protein import classify as cl
from dnainvent.protein.orf import best_orf, find_orfs
from i18n import LOCALES_DIR, available_languages

AMINO = "ACDEFGHIKLMNPQRSTVWY"


# ---------- Utilidades: puntuar un alineamiento ya hecho y enumerar todos ----------

def score_of(first: str, second: str, scoring) -> int:
    """Puntuación de un alineamiento dado, calculada columna a columna (sin tablas)."""
    total, previous = 0, None          # previous: "X" (hueco abajo), "Y" (hueco arriba) o None
    for a, b in zip(first, second):
        if a == alignment.GAP:
            total += scoring.gap_extend if previous == "Y" else scoring.gap_open
            previous = "Y"
        elif b == alignment.GAP:
            total += scoring.gap_extend if previous == "X" else scoring.gap_open
            previous = "X"
        else:
            total += scoring.pair(a, b)
            previous = None
    return total


def all_alignments(a: str, b: str):
    """Todas las formas de alinear a y b (solo para secuencias muy cortas)."""
    @lru_cache(maxsize=None)
    def go(i, j):
        if i == len(a) and j == len(b):
            return [("", "")]
        out = []
        if i < len(a) and j < len(b):
            out += [(a[i] + x, b[j] + y) for x, y in go(i + 1, j + 1)]
        if i < len(a):
            out += [(a[i] + x, "-" + y) for x, y in go(i + 1, j)]
        if j < len(b):
            out += [("-" + x, b[j] + y) for x, y in go(i, j + 1)]
        return out
    return go(0, 0)


# ---------- BLOSUM62 ----------

def test_blosum62_is_symmetric_and_complete():
    for x in BLOSUM62_ORDER:
        for y in BLOSUM62_ORDER:
            assert BLOSUM62[x][y] == BLOSUM62[y][x], (x, y)


def test_blosum62_known_values():
    assert BLOSUM62["W"]["W"] == 11      # el triptófano es el más raro: coincidir vale mucho
    assert BLOSUM62["C"]["C"] == 9
    assert BLOSUM62["I"]["V"] == 3       # parecidos (pequeños e hidrófobos)
    assert BLOSUM62["W"]["D"] == -4      # muy distintos


# ---------- Alineamiento de proteínas ----------

def test_protein_identical():
    result = alignment.align("MKTAYIAKQR", "MKTAYIAKQR", kind="protein")
    assert result.identity == 100
    assert result.kind == "protein"
    assert result.score == sum(BLOSUM62[c][c] for c in "MKTAYIAKQR")


def test_protein_similar_substitution_scores_positive():
    # M-M (5) + I-V (3) + V-V (4): cambiar I por V apenas penaliza.
    assert alignment.align("MIV", "MVV", kind="protein").score == 12


def test_affine_gaps_make_one_long_gap():
    # Con huecos afines, los 5 aminoácidos que faltan forman un solo hueco, no varios.
    result = alignment.align("MKTAYIAKQRQISFVKSH", "MKTAYIQISFVKSH", kind="protein")
    assert result.second.count("-") == 4
    assert "----" in result.second


def test_protein_rejects_invalid_letters():
    with pytest.raises(InvalidSequenceError) as info:
        alignment.align("MK1", "MK", kind="protein")
    assert info.value.code == "invalid_characters"
    assert info.value.params["kind"] == "protein"


def test_unknown_kind_is_rejected():
    with pytest.raises(ValueError):
        alignment.align("ACGT", "ACGT", kind="rna")


def test_dna_is_still_the_default():
    assert alignment.align("ACGT", "AGGT").kind == "dna"
    with pytest.raises(InvalidSequenceError):
        alignment.align("MKV", "MKV")          # sin kind="protein" se valida como ADN


@pytest.mark.parametrize("kind", ["dna", "protein"])
def test_alignment_is_optimal_by_brute_force(kind):
    """Para secuencias cortas se prueban TODOS los alineamientos: ninguno puntúa más."""
    rng = random.Random(7)
    letters = "ACGT" if kind == "dna" else AMINO
    scoring = alignment.KINDS[kind][1]
    for _ in range(40):
        a = "".join(rng.choice(letters) for _ in range(rng.randint(1, 5)))
        b = "".join(rng.choice(letters) for _ in range(rng.randint(1, 5)))
        best = max(score_of(x, y, scoring) for x, y in all_alignments(a, b))
        assert alignment.align(a, b, kind=kind).score == best, (a, b)


@pytest.mark.parametrize("kind", ["dna", "protein"])
def test_returned_alignment_has_the_reported_score(kind):
    """El camino reconstruido suma exactamente la puntuación que se devuelve."""
    rng = random.Random(3)
    letters = "ACGT" if kind == "dna" else AMINO
    scoring = alignment.KINDS[kind][1]
    for _ in range(50):
        a = "".join(rng.choice(letters) for _ in range(rng.randint(1, 60)))
        b = "".join(rng.choice(letters) for _ in range(rng.randint(1, 60)))
        r = alignment.align(a, b, kind=kind)
        assert r.first.replace("-", "") == a and r.second.replace("-", "") == b
        assert score_of(r.first, r.second, scoring) == r.score


def test_custom_gap_costs():
    cheap = alignment.MatrixScoring(gap_open=-1, gap_extend=-1)
    result = alignment.align("WWWW", "DDDD", cheap, kind="protein")
    assert result.mismatches == 0              # con huecos baratos, evita W frente a D


def test_pairs_maps_positions():
    result = alignment.align("ACGTACGTAC", "TACGTACGTAC")
    assert result.pairs()[:3] == [(0, 1), (1, 2), (2, 3)]
    assert len(result.pairs()) == result.length - result.gaps


# ---------- ORF ----------

def test_simple_forward_orf():
    orf = best_orf("ATGAAATAG")
    assert (orf.protein, orf.stop_codon, orf.start, orf.end) == ("MK", "TAG", 0, 9)
    assert orf.label == "+1" and orf.is_complete


def test_orf_matches_translate_on_forward_strand():
    dna = "GC" + "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACTGA" + "GG"
    orf = best_orf(dna)
    assert orf.protein == sequences.translate(sequences.transcribe(dna)).protein
    assert orf.protein == "MVHLTPEEKSAVTALWGKVN"   # inicio de la β-globina humana
    assert orf.label == "+3"


def test_reverse_strand_orf_and_coordinates():
    gene = "ATGAAACCCGGGTAA"
    dna = "TT" + sequences.reverse_complement(gene) + "C"
    orf = best_orf(dna)
    assert orf.strand == -1 and orf.protein == "MKPG"
    assert sequences.reverse_complement(dna[orf.start:orf.end]) == gene


def test_incomplete_orf_runs_to_the_end():
    orf = best_orf("ATGAAACC")             # sin parada y con un codón a medias
    assert orf.protein == "MK"
    assert orf.stop_codon is None and not orf.is_complete
    assert orf.end == 6


def test_inner_atg_does_not_start_a_new_orf():
    orfs = [o for o in find_orfs("ATGATGAAATAA") if o.strand == 1 and o.frame == 0]
    assert [o.protein for o in orfs] == ["MMK"]


def test_no_start_codon():
    assert find_orfs("CCCCCCCC") == []
    with pytest.raises(InvalidSequenceError) as info:
        best_orf("CCCCCCCC")
    assert info.value.code == "no_orf"


def test_min_length():
    dna = "ATGTAA" + "ATGAAACCCGGGTTTTAA"
    assert [o.protein for o in find_orfs(dna, min_length=3)] == ["MKPGF"]
    with pytest.raises(InvalidParameterError) as info:
        find_orfs(dna, min_length=0)
    assert info.value.code == "min_length_too_small"


def test_complete_orf_wins_a_tie():
    # ATG CCC TAA -> "MP" completo; ATG AAA -> "MK" sin parada. La otra hebra no tiene ATG.
    orfs = find_orfs("ATGCCCTAAATGAAA")
    assert [(o.protein, o.is_complete) for o in orfs] == [("MP", True), ("MK", False)]


def test_every_orf_is_consistent_with_its_coordinates():
    """En secuencias al azar, cada ORF se puede releer desde sus coordenadas."""
    dna_codons = {k.replace("U", "T"): v for k, v in CODON_TABLE.items()}
    rng = random.Random(11)
    for _ in range(30):
        dna = "".join(rng.choice("ACGT") for _ in range(rng.randint(3, 300)))
        orfs = find_orfs(dna)
        assert [o.length for o in orfs] == sorted((o.length for o in orfs), reverse=True)
        for orf in orfs:
            piece = dna[orf.start:orf.end]
            if orf.strand < 0:
                piece = sequences.reverse_complement(piece)
            codons = [piece[i:i + 3] for i in range(0, len(piece), 3)]
            amino = "".join(dna_codons[c] for c in codons)
            assert piece.startswith("ATG")
            if orf.is_complete:
                assert amino == orf.protein + "_"
            else:
                assert amino == orf.protein and "_" not in amino
            # la pauta es la posición de inicio módulo 3 dentro de su hebra
            first = orf.start if orf.strand > 0 else len(dna) - orf.end
            assert first % 3 == orf.frame


# ---------- Clasificación ----------

def hit(identity, query=100.0, subject=100.0, evalue=1e-50, accession="P00001"):
    return cl.Hit(accession, identity, query, subject, evalue)


@pytest.mark.parametrize("h, expected", [
    (hit(100), cl.Category.IDENTICAL),
    (hit(99.5), cl.Category.VERY_SIMILAR),
    (hit(90), cl.Category.VERY_SIMILAR),
    (hit(89.9), cl.Category.DISTANT_RELATIVE),
    (hit(30), cl.Category.DISTANT_RELATIVE),
    (hit(25), cl.Category.INVENTED),                    # zona crepuscular
    (hit(100, subject=40), cl.Category.VERY_SIMILAR),   # trozo de una proteína real
    (hit(100, query=50), cl.Category.INVENTED),         # cubre poco de la consulta
    (hit(100, evalue=0.5), cl.Category.INVENTED),       # podría ser casualidad
])
def test_single_hit(h, expected):
    assert cl.classify([h]).category == expected


def test_no_hits_is_invented():
    result = cl.classify([])
    assert result.category == cl.Category.INVENTED
    assert result.hit is None and not result.is_known


def test_best_hit_wins():
    hits = [hit(45, accession="DISTANT"), hit(100, accession="EXACT"),
            hit(95, accession="VARIANT")]
    result = cl.classify(hits)
    assert result.category == cl.Category.IDENTICAL
    assert result.hit.accession == "EXACT"


def test_tie_breaks_by_identity_then_evalue():
    hits = [hit(92, evalue=1e-30, accession="A"), hit(96, evalue=1e-20, accession="B"),
            hit(96, evalue=1e-40, accession="C")]
    assert cl.classify(hits).hit.accession == "C"


def test_preferred_hit_breaks_exact_ties():
    hits = [hit(100, accession="CHIMP"), hit(100, accession="HUMAN"), hit(100, accession="BONOBO")]
    assert cl.classify(hits).hit.accession == "CHIMP"              # sin preferencia: el 1.º
    hits[1] = cl.Hit("HUMAN", 100, 100, 100, 1e-50, preferred=True)
    assert cl.classify(hits).hit.accession == "HUMAN"
    better = [hit(100, accession="CHIMP", evalue=1e-60), hits[1]]
    assert cl.classify(better).hit.accession == "CHIMP"            # solo desempata


def test_custom_thresholds():
    strict = cl.Thresholds(very_similar_identity=98, distant_identity=50)
    assert cl.classify([hit(95)], strict).category == cl.Category.DISTANT_RELATIVE
    assert cl.classify([hit(45)], strict).category == cl.Category.INVENTED


@pytest.mark.parametrize("bad", [
    dict(identity=101), dict(query_coverage=-1), dict(subject_coverage=200), dict(evalue=-1),
])
def test_hit_rejects_impossible_values(bad):
    values = dict(accession="X", identity=50, query_coverage=100, subject_coverage=100,
                  evalue=1e-10) | bad
    with pytest.raises(ValueError):
        cl.Hit(**values)


def test_thresholds_must_be_ordered():
    with pytest.raises(ValueError):
        cl.Thresholds(very_similar_identity=30, distant_identity=90)


@pytest.mark.parametrize("language", sorted(available_languages()))
def test_every_category_has_texts(language):
    data = json.loads((LOCALES_DIR / f"{language}.json").read_text(encoding="utf-8"))
    for category in cl.Category:
        node = data
        for part in category.i18n_key.split("."):
            node = node[part]
        assert node["name"] and node["description"]
