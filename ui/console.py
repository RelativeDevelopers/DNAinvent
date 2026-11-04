"""Interfaz de consola. Aquí (y solo aquí) hay input(), print() y colores.

Ningún texto está escrito en este archivo: todo sale de i18n/locales/<idioma>.json
a través de tr("clave"). Para añadir un idioma basta con crear otro JSON.
"""

import re
import sys

from colorama import Fore, Style, init

from dnainvent import DNAinventError, alignment, bits, crispr, sequences, tools
from dnainvent.constants import AMINO_ACID_CODES
from dnainvent.database import default_path, open_database, queries
from dnainvent.database.clinvar import ClinVarProvider
from dnainvent.database.ncbi import HUMAN, NCBIProvider
from dnainvent.database.variants import locate, side_by_side, three_letter
from dnainvent.protein import confidence, pipeline
from dnainvent.protein.cache import key_for
from dnainvent.protein.classify import Category, Thresholds
from dnainvent.protein.clients.afdb import AlphaFoldDB
from dnainvent.protein.clients.blast import BlastClient, blast_email
from dnainvent.protein.compare import moved
from dnainvent.protein.predictors.esmfold import ESMFold
from dnainvent.validation import validate_bits, validate_dna, validate_rna
from i18n import DEFAULT_LANGUAGE, Translator, available_languages

from . import launcher, settings

BASE_COLORS = {"A": Fore.YELLOW, "C": Fore.BLUE, "G": Fore.GREEN, "T": Fore.RED}
CATEGORY_COLORS = {Category.IDENTICAL: Fore.GREEN, Category.VERY_SIMILAR: Fore.CYAN,
                   Category.DISTANT_RELATIVE: Fore.YELLOW, Category.INVENTED: Fore.MAGENTA}
OTHER_HITS = 5            # parecidos que se enseñan además del principal
MAX_MUTATIONS = 15        # mutaciones que se enumeran (el resto se cuentan)
MOST_DIFFERENT = 5        # zonas de forma más distinta que se enseñan
DEVIATION_SHOWN = 2.0     # Å: por debajo de esto no se considera que la zona cambie
PERMUTATION_WARNING = 10_000
ALIGNMENT_WARNING = 4_000_000   # unas 2.000 x 2.000 bases
LINE_WIDTH = 60
VARIANTS_SHOWN = 20       # variantes que se enumeran al buscar
DISEASES_SHOWN = 10       # enfermedades en la ficha de un gen

# Clasificación de ClinVar -> clave de i18n (database.significance.*)
SIGNIFICANCE_KEYS = {
    "pathogenic": "pathogenic", "likely pathogenic": "likely_pathogenic",
    "pathogenic/likely pathogenic": "pathogenic_likely_pathogenic",
    "uncertain significance": "uncertain", "likely benign": "likely_benign",
    "benign": "benign", "benign/likely benign": "benign_likely_benign",
    "conflicting classifications of pathogenicity": "conflicting",
    "protective": "protective", "risk factor": "risk_factor", "other": "other",
    "drug response": "drug_response", "association": "association", "affects": "affects",
    "not provided": "not_provided", "oncogenic": "oncogenic", "likely oncogenic": "likely_oncogenic",
    "pathogenic, low penetrance": "pathogenic_low_penetrance",
    "likely pathogenic, low penetrance": "likely_pathogenic_low_penetrance",
    "established risk allele": "established_risk_allele", "likely risk allele": "likely_risk_allele",
    "uncertain risk allele": "uncertain_risk_allele",
    "tier i - strong": "tier_1", "tier ii - potential": "tier_2", "tier iii - unknown": "tier_3",
    "tier iv - benign/likely benign": "tier_4",
    "no classification for the single variant": "no_single_classification",
}
SIGNIFICANCE_COLORS = {"pathogenic": Fore.RED, "likely_pathogenic": Fore.RED,
                       "pathogenic_likely_pathogenic": Fore.RED, "oncogenic": Fore.RED,
                       "likely_oncogenic": Fore.RED, "pathogenic_low_penetrance": Fore.RED,
                       "likely_pathogenic_low_penetrance": Fore.RED, "tier_1": Fore.RED, "uncertain": Fore.YELLOW,
                       "conflicting": Fore.YELLOW, "likely_benign": Fore.GREEN,
                       "benign": Fore.GREEN, "benign_likely_benign": Fore.GREEN,
                       "protective": Fore.CYAN, "risk_factor": Fore.MAGENTA}

tr = Translator(DEFAULT_LANGUAGE)


# ---------- Utilidades de entrada y salida ----------

def error(text: str):
    print(Fore.RED + text + Style.RESET_ALL)


def warning(text: str):
    print(Fore.YELLOW + text + Style.RESET_ALL)


def ask_sequence(prompt_key: str, validator, hint_key: str | None = None) -> str:
    """Pide una secuencia hasta que sea válida."""
    while True:
        text = input(tr(prompt_key))
        try:
            seq = validator(text)
            print(tr("common.sequence_length", length=len(seq)))
            return seq
        except DNAinventError as exc:
            message = tr.error(exc)
            if hint_key and "T" in text.upper():
                message += " " + tr(hint_key)
            error(message)


def ask_dna(prompt_key: str = "prompts.dna") -> str:
    return ask_sequence(prompt_key, validate_dna)


def ask_rna() -> str:
    # Si el usuario escribe ADN (con T), se le sugiere transcribirlo antes.
    return ask_sequence("prompts.rna", validate_rna, hint_key="translate.is_dna_hint")


def ask_int(prompt: str, minimum: int = 0, maximum: int | None = None,
            default: int | None = None) -> int:
    while True:
        text = input(prompt).strip()
        if not text and default is not None:
            return default
        try:
            value = int(text)
            if value >= minimum and (maximum is None or value <= maximum):
                return value
        except ValueError:
            pass
        if maximum is None:
            error(tr("common.int_minimum", minimum=minimum))
        else:
            error(tr("common.int_between", minimum=minimum, maximum=maximum))


def confirm(question_key: str, **params) -> bool:
    answer = input(tr(question_key, **params) + tr("common.confirm_suffix")).strip().lower()
    return answer in tr("common.yes_answers").split(",")


def print_labeled(rows: list[tuple[str, str]], indent: str = ""):
    """Imprime pares (clave de etiqueta, valor) con las etiquetas alineadas."""
    labels = [tr(key) for key, _ in rows]
    width = max(len(label) for label in labels) + 1
    for label, (_, value) in zip(labels, rows):
        print(f"{indent}{label:<{width}}{value}")


def highlight(seq: str, positions: list[int], length: int) -> str:
    """Resalta en cian los fragmentos que empiezan en positions."""
    out, last = [], 0
    for pos in positions:
        out.append(seq[last:pos])
        out.append(Fore.CYAN + seq[pos:pos + length] + Style.RESET_ALL)
        last = pos + length
    out.append(seq[last:])
    return "".join(out)


def color_dna(dna: str) -> str:
    return "".join(BASE_COLORS[b] + b for b in dna) + Style.RESET_ALL


def color_bits(bit_string: str) -> str:
    """Bits agrupados de 2 en 2 con el color de la base que representan."""
    return " ".join(BASE_COLORS[base] + block
                    for block, base in zip(bits.split_blocks(bit_string),
                                           sequences.from_binary(bit_string))) + Style.RESET_ALL


# ---------- Menú principal ----------

def option_complement():
    dna = ask_dna()
    print("\n" + tr("complement.helix"))
    print(f"  5'-{dna}-3'")
    print(f"     {'|' * len(dna)}")
    print(f"  3'-{sequences.complement(dna)}-5'")
    print("\n" + tr("complement.reverse"))
    print(Fore.CYAN + f"  5'-{sequences.reverse_complement(dna)}-3'" + Style.RESET_ALL)
    print()
    if confirm("microscope.offer"):
        explore_in_3d(dna)


def option_transcribe():
    print(sequences.transcribe(ask_dna()))


def option_translate():
    rna = ask_rna()
    result = sequences.translate(rna)

    if not result.has_start:
        warning(tr("translate.no_start"))
        return

    # Marco de lectura en cian y codón de parada en rojo.
    stop_len = 3 if result.is_complete else 0
    orf_end = result.end - stop_len
    print(rna[:result.start]
          + Fore.CYAN + rna[result.start:orf_end]
          + Fore.RED + rna[orf_end:result.end] + Style.RESET_ALL
          + rna[result.end:])

    print("\n" + tr("translate.start_at", position=result.start + 1))
    if result.is_complete:
        print(tr("translate.stop_at", codon=result.stop_codon, position=orf_end + 1))
    else:
        warning(tr("translate.no_stop"))

    print("\n" + tr("translate.protein", count=len(result.protein), protein=result.protein))
    for aa, count in sequences.count_amino_acids(result.protein).most_common():
        print(f"  {tr.amino_acid(aa):<22} x{count}")


def paint_target(dna: str, target: crispr.Target) -> str:
    """Protoespaciador en cian, PAM en magenta y el corte marcado con una barra roja."""
    out = []
    for index, base in enumerate(dna):
        if index == target.cut:
            out.append(Fore.RED + Style.BRIGHT + "|" + Style.RESET_ALL)
        if target.start <= index < target.end:
            out.append(Fore.CYAN + base)
        elif target.pam_start <= index < target.pam_end:
            out.append(Fore.MAGENTA + base)
        else:
            out.append(Style.RESET_ALL + base)
    return "".join(out) + Style.RESET_ALL


def show_protein_effect(before: str, after: str):
    """Traduce la hebra antes y después de editar para ver el efecto en la proteína."""
    old = sequences.translate(sequences.transcribe(before))
    new = sequences.translate(sequences.transcribe(after))
    if not old.has_start and not new.has_start:
        print("  " + tr("crispr.no_atg"))
        return
    width = max(len(tr("crispr.before")), len(tr("crispr.after"))) + 1
    print(f"  {tr('crispr.before'):<{width}}{old.protein or '-'}")
    print(f"  {tr('crispr.after'):<{width}}{new.protein or '-'}")
    if old.protein == new.protein:
        print("  " + tr("crispr.unchanged"))


def describe_edit(edit: crispr.Edit) -> str:
    # Posiciones contadas desde 1, como en biología.
    first, last = edit.position + 1, edit.position + len(edit.deleted)
    if edit.kind == "insertion":
        return tr("crispr.insertion", count=len(edit.inserted), bases=edit.inserted,
                  before=edit.position, after=edit.position + 1)
    if edit.kind == "deletion":
        return tr("crispr.deletion", count=len(edit.deleted), bases=edit.deleted,
                  first=first, last=last)
    return tr("crispr.substitution", deleted=edit.deleted, inserted=edit.inserted,
              first=first, last=last)


def option_crispr():
    dna = ask_dna()
    guide_length = ask_int(tr("prompts.guide_length", default=crispr.GUIDE_LENGTH),
                           minimum=crispr.CUT_OFFSET + 1, default=crispr.GUIDE_LENGTH)

    # 1) Buscar dianas
    targets = crispr.find_targets(dna, guide_length)
    if not targets:
        warning(tr("crispr.no_targets", length=guide_length))
        return

    print("\n" + tr("crispr.found", count=len(targets)))
    for number, t in enumerate(targets, start=1):
        print(f"  {number:>2}] " + tr("crispr.target", strand=t.strand, cut=f"{t.cut:>4}",
                                      guide=t.guide_rna, pam=t.pam))
    choice = ask_int(tr("prompts.choose_target"), minimum=1, maximum=len(targets))
    target = targets[choice - 1]

    # 2) Cortar
    print("\n" + tr("crispr.cut_here"))
    print("  " + paint_target(dna, target))

    # 3) Reparar
    print("\n" + tr("crispr.repair_question"))
    print("  1] " + tr("crispr.nhej"))
    print("  2] " + tr("crispr.hdr"))
    repair = ask_int(tr("common.choose_option"), minimum=1, maximum=2)

    if repair == 1:
        edit = crispr.repair_nhej(dna, target.cut)
    else:
        insert = ask_dna("prompts.insert")
        edit = crispr.repair_hdr(dna, target.cut, insert)
        left, _, right = crispr.donor_template(dna, target.cut, insert)
        print("\n" + tr("crispr.donor"))
        print("  " + Style.DIM + left + Style.RESET_ALL + Fore.GREEN + insert
              + Style.RESET_ALL + Style.DIM + right + Style.RESET_ALL)

    # Resultado
    print("\n" + tr("crispr.result", change=describe_edit(edit)))
    show_alignment(alignment.align(edit.original, edit.edited))
    if edit.frameshift:
        warning(tr("crispr.frameshift"))
    else:
        print(tr("crispr.in_frame"))
    print("\n" + tr("crispr.protein_effect"))
    show_protein_effect(edit.original, edit.edited)


def option_binary():
    dna = ask_dna()
    bit_string = sequences.to_binary(dna)
    print(color_bits(bit_string))
    print("\n" + tr("binary.summary", bases=len(dna), bits=len(bit_string),
                    bytes=tr.decimal(len(bit_string) / 8, 3, trim=True)))
    print(tr("binary.text_size", bits=len(dna) * 8))


def option_from_binary():
    bit_string = ask_sequence("prompts.binary", validate_bits)
    print(color_dna(sequences.from_binary(bit_string)))


def option_microscope():
    explore_in_3d(ask_dna())


def explore_in_3d(dna: str):
    """Elige el tramo (si la secuencia es larga) y abre el microscopio en otra ventana."""
    from dnainvent.scene.limits import RECOMMENDED_PAIRS   # no necesita NumPy

    missing = launcher.missing_dependencies()
    if missing:
        error(tr("microscope.missing", packages=missing))
        print(tr("microscope.install"))
        return
    start, count = 0, len(dna)
    if len(dna) > RECOMMENDED_PAIRS:
        print(tr("microscope.long", length=len(dna), maximum=RECOMMENDED_PAIRS))
        start = ask_int(tr("microscope.ask_start"), minimum=1, maximum=len(dna), default=1) - 1
        count = min(RECOMMENDED_PAIRS, len(dna) - start)
    print(tr("microscope.window", first=start + 1, last=start + count, length=len(dna)))
    launcher.open_microscope(dna, start, count, tr.language)
    print(tr("microscope.opening"))


def option_protein():
    """Del ADN a la proteína: marco de lectura, BLAST, caso, estructura natural (AlphaFold DB)
    y, si no es idéntica a una real, la predicción de la tuya (ESMFold) y la comparación."""
    dna = ask_dna()
    email = blast_email()
    print(tr("protein.email.using", email=email))
    try:
        report = pipeline.analyze(dna, BlastClient(email), AlphaFoldDB(), ESMFold(),
                                  progress=ProgressPrinter())
    except DNAinventError as exc:
        print()
        error(tr.error(exc))
        return
    print()
    show_protein_report(report)
    if report.structure or report.prediction:
        if confirm("protein.save_offer"):
            save_structures(report)
        if confirm("protein.microscope_offer"):
            explore_protein_in_3d(report)


def explore_protein_in_3d(report: pipeline.ProteinReport):
    """Prepara la sesión (la natural y/o la tuya) y abre el microscopio en otra ventana."""
    from dnainvent.protein.session import panels_for, write_session
    from dnainvent.scene.limits import RECOMMENDED_RESIDUES

    missing = launcher.missing_dependencies()
    if missing:
        error(tr("microscope.missing", packages=missing))
        print(tr("microscope.install"))
        return
    main = report.structure if report.structure is not None else report.prediction
    length = len(main.residues())
    start = 0
    if length > RECOMMENDED_RESIDUES:
        print(tr("protein.microscope.long", length=length, maximum=RECOMMENDED_RESIDUES))
        start = ask_int(tr("protein.microscope.ask_start"), minimum=1, maximum=length,
                        default=1) - 1
    count = min(RECOMMENDED_RESIDUES, length - start)
    natural = (f"{report.best.description} · {report.best.accession}" if report.best
               else report.entry.model_id if report.entry else "")
    labels = {"natural": natural,
              "yours": tr("protein.microscope.yours", predictor=report.predictor)}
    session = write_session(panels_for(report, labels, start, count))
    print(tr("protein.microscope.window", first=start + 1, last=start + count, length=length))
    launcher.open_protein_microscope(session, tr.language)
    print(tr("microscope.opening"))


class ProgressPrinter:
    """Cada paso en su línea; mientras BLAST trabaja, un punto por cada consulta."""

    def __init__(self):
        self.dots = False

    def __call__(self, key: str, data: dict):
        if key == "protein.stage.waiting":
            print(("." if self.dots else " ."), end="", flush=True)
            self.dots = True
            return
        self.dots = False
        if "database" in data:
            data = {**data, "database": tr(f"protein.databases.{data['database']}")}
        print("\n" + tr(key, **data), end="", flush=True)


def show_protein_report(report: pipeline.ProteinReport):
    orf = report.orf
    stop = (tr("protein.report.complete", codon=orf.stop_codon) if orf.is_complete
            else tr("protein.report.incomplete"))
    print("\n" + tr("protein.report.orf", frame=orf.label, start=orf.start + 1, end=orf.end,
                    length=report.dna_length, count=orf.length, stop=stop))
    print(tr("protein.report.protein"))
    for start in range(0, len(report.protein), LINE_WIDTH):
        print("  " + report.protein[start:start + LINE_WIDTH])

    color = CATEGORY_COLORS[report.category]
    key = report.category.i18n_key
    limits = Thresholds()
    thresholds = {"high": tr.decimal(limits.very_similar_identity, 0),
                  "low": tr.decimal(limits.distant_identity, 0)}
    print("\n" + tr("protein.report.category") + " " + color + Style.BRIGHT
          + tr(f"{key}.name") + Style.RESET_ALL)
    print("  " + tr(f"{key}.description", **thresholds))

    if report.best is not None:
        print("\n" + tr("protein.report.best"))
        print("  " + describe_hit(report.best))
        print("  " + tr("protein.report.identity", identity=tr.decimal(report.best.identity, 1),
                        query=tr.decimal(report.best.query_coverage, 0),
                        subject=tr.decimal(report.best.subject_coverage, 0),
                        evalue=tr.sci(report.best.evalue)))
        if report.category is not Category.IDENTICAL:
            show_blast_alignment(report.best)
        others = [h for h in report.hits if h is not report.best][:OTHER_HITS]
        if others:
            print("\n" + tr("protein.report.others"))
            for hit in others:
                print(f"  {tr.decimal(hit.identity, 1):>5} %  {describe_hit(hit)}")
    elif report.database is not None:
        print("\n" + tr("protein.report.none",
                        database=tr(f"protein.databases.{report.database}")))

    if report.structure is not None:
        show_structure_summary(report)
    if report.prediction is not None:
        show_prediction_summary(report)
    if report.comparison is not None:
        show_comparison(report)
    if report.warnings:
        print("\n" + tr("protein.report.warnings"))
        for warning_key, data in report.warnings:
            if "plddt" in data:
                data = {**data, "plddt": tr.decimal(data["plddt"], 1)}
            warning("  - " + tr(warning_key, **data))


def describe_hit(hit) -> str:
    parts = [Fore.CYAN + hit.accession + Style.RESET_ALL, hit.description, hit.organism]
    if hit.gene:
        parts.append(tr("protein.report.gene", gene=hit.gene))
    return " · ".join(part for part in parts if part)


def show_blast_alignment(hit):
    print("\n" + tr("protein.report.alignment"))
    for start in range(0, len(hit.query_aligned), LINE_WIDTH):
        end = start + LINE_WIDTH
        print("  " + hit.query_aligned[start:end])
        print("  " + hit.markup[start:end])
        print("  " + hit.hit_aligned[start:end])
        print()


def show_confidence(structure):
    for _, band, _ in confidence.BANDS:
        print(f"  {tr(f'protein.confidence.{band}')}: {confidence.counts(structure)[band]}")


def show_structure_summary(report: pipeline.ProteinReport):
    entry = report.entry
    print("\n" + tr("protein.report.structure", model=entry.model_id, version=entry.version,
                    plddt=tr.decimal(confidence.mean(report.structure), 1)))
    show_confidence(report.structure)
    print(tr("protein.report.citation"))


def show_prediction_summary(report: pipeline.ProteinReport):
    print("\n" + tr("protein.report.prediction", predictor=report.predictor,
                    plddt=tr.decimal(confidence.mean(report.prediction), 1)))
    show_confidence(report.prediction)
    print(tr("protein.report.prediction_citation", predictor=report.predictor))


def show_comparison(report: pipeline.ProteinReport):
    c = report.comparison
    print("\n" + tr("protein.report.comparison", accession=report.best.accession))
    if c.differences:
        labels = [d.label for d in c.differences[:MAX_MUTATIONS]]
        if len(c.differences) > MAX_MUTATIONS:
            labels.append(tr("protein.report.differences_more",
                             count=len(c.differences) - MAX_MUTATIONS))
        key = ("protein.report.differences_distant"
               if report.category is Category.DISTANT_RELATIVE else "protein.report.differences")
        print("  " + tr(key, count=len(c.differences), list=", ".join(labels)))
    else:
        print("  " + tr("protein.report.no_differences"))
    if c.insertions or c.deletions:
        print("  " + tr("protein.report.indels", insertions=c.insertions, deletions=c.deletions))
    if report.category is Category.DISTANT_RELATIVE:
        print("  " + tr("protein.report.conserved", count=len(c.conserved), total=len(c.pairs),
                        percent=tr.decimal(c.identity, 1)))
    if c.rmsd is not None:
        meaning = ("rmsd_same" if c.rmsd < 1 else "rmsd_close" if c.rmsd <= 2
                   else "rmsd_different")
        print("  " + tr("protein.report.rmsd", rmsd=tr.decimal(c.rmsd, 2),
                        count=c.superposition.count, meaning=tr(f"protein.report.{meaning}")))
        natural = report.structure.residues()
        worst = sorted((d, j) for j, d in c.deviations.items() if d >= DEVIATION_SHOWN)
        worst = sorted(worst[-MOST_DIFFERENT:], key=lambda item: item[1])
        if worst:
            labels = [f"{AMINO_ACID_CODES.get(natural[j].letter, 'Xaa')}{j + 1} "
                      f"({tr.decimal(d, 1)} Å)" for d, j in worst]
            print("  " + tr("protein.report.most_different", list=", ".join(labels)))
    if c.differences:
        print("  " + tr("protein.report.numbering"))


def save_structures(report: pipeline.ProteinReport):
    """La natural y la tuya (superpuesta sobre la natural), con hidrógenos."""
    saved = []
    if report.structure is not None:
        saved.append(save_with_hydrogens(report.structure,
                                         f"{report.entry.model_id}_DNAinvent.pdb"))
    if report.prediction is not None:
        prediction = report.prediction
        if report.comparison is not None and report.comparison.superposition is not None:
            prediction = moved(prediction, report.comparison.superposition)
        name = f"{report.predictor}_{key_for(report.protein)[:8]}_DNAinvent.pdb"
        saved.append(save_with_hydrogens(prediction, name))
    if len(saved) == 2:
        print(tr("protein.saved_pair"))


def save_with_hydrogens(structure, filename: str):
    from pathlib import Path

    from dnainvent.protein.hydrogens import add_hydrogens
    from dnainvent.protein.pdb import write_pdb

    result = add_hydrogens(structure.residues())
    atoms = [atom for residue in result.residues for atom in residue.atoms.values()]
    path = Path.cwd() / filename
    path.write_text(write_pdb(atoms, structure.title + " + H (DNAinvent)"), encoding="utf-8")
    print(tr("protein.report.atoms", heavy=len(atoms) - result.added, hydrogens=result.added))
    print(tr("protein.saved", path=path))
    return path


# ---------- Base de datos de genes reales ----------

_database = None


def database():
    """La base de datos se abre la primera vez que hace falta y queda abierta."""
    global _database
    if _database is None:
        _database = open_database()
    return _database


class DatabaseProgress:
    """Cada paso en su línea; los que van por lotes (200 de 1.974...) se reescriben en la misma."""

    def __init__(self):
        self.same_line = None

    def __call__(self, key: str, data: dict):
        text = tr(key, **data)
        if "done" in data:
            prefix = "\r" if self.same_line == key else ("\n" if self.same_line else "")
            print(prefix + "  " + text, end="", flush=True)
            self.same_line = key
        else:
            print(("\n" if self.same_line else "") + "  " + text, flush=True)
            self.same_line = None

    def finish(self):
        if self.same_line:
            print()
        self.same_line = None


def significance(value: str | None) -> str:
    """Clasificación de ClinVar en el idioma del usuario (y en color)."""
    if not value:
        return tr("database.significance.none")

    def one(part: str) -> str:
        key = SIGNIFICANCE_KEYS.get(part)
        if key is None:
            return part
        return SIGNIFICANCE_COLORS.get(key, "") + tr(f"database.significance.{key}") + Style.RESET_ALL

    if value in SIGNIFICANCE_KEYS:
        return one(value)
    # ClinVar a veces junta varias: "pathogenic; other", "pathogenic/likely pathogenic/
    # pathogenic, low penetrance". Se traduce cada trozo y se dejan los separadores.
    return "".join(one(part) if i % 2 == 0 else part
                   for i, part in enumerate(re.split(r"(; |/)", value)))


def ask_gene():
    """Pide un símbolo hasta que sea un gen guardado (o se deje vacío). Devuelve la fila o None."""
    if not queries.genes(database()):
        warning(tr("database.empty"))
        return None
    while True:
        symbol = input(tr("prompts.gene_symbol")).strip()
        if not symbol:
            return None
        gene = queries.find_gene(database(), symbol)
        if gene is not None:
            return gene
        warning(tr("database.not_saved", symbol=symbol.upper()))


def db_import():
    symbol = input(tr("prompts.gene_symbol")).strip()
    if not symbol:
        return
    taxon = ask_int(tr("prompts.taxon"), minimum=1, default=HUMAN)
    conn, progress = database(), DatabaseProgress()
    try:
        gene_id = NCBIProvider().import_gene(conn, symbol, taxon, progress)
        gene = conn.execute("SELECT * FROM genes WHERE id = ?", (gene_id,)).fetchone()
        print(Fore.GREEN + tr("database.imported_gene", symbol=gene["symbol"], name=gene["name"])
              + Style.RESET_ALL)
        if taxon != HUMAN:
            print(tr("database.clinvar_human_only"))
            return
        if not confirm("database.clinvar_offer"):
            return
        stats = ClinVarProvider().import_variants(conn, gene["symbol"], progress=progress)
    except DNAinventError as exc:
        progress.finish()
        error(tr.error(exc))
        return
    progress.finish()
    print()
    print_labeled([
        ("database.stats.found", tr.number(stats.found)),
        ("database.stats.saved", tr("database.stats.saved_value", count=stats.saved,
                                    pathogenic=stats.pathogenic)),
        ("database.stats.with_conditions", tr.number(stats.with_conditions)),
        ("database.stats.not_on_mrna", tr.number(stats.not_on_mrna)),
        ("database.stats.other_gene", tr.number(stats.other_gene)),
        ("database.stats.new_diseases", tr.number(stats.new_diseases)),
    ], indent="  ")
    if stats.saved and stats.with_conditions < stats.saved / 2:
        warning(tr("database.stats.few_conditions", count=stats.with_conditions, saved=stats.saved))
    if stats.condition_batches_failed:
        warning(tr("database.stats.batches_failed", count=stats.condition_batches_failed))
    if stats.unexpected_records:
        warning(tr("database.stats.unexpected", count=stats.unexpected_records))
    print(Style.DIM + tr("database.citation") + Style.RESET_ALL)


def db_gene():
    gene = ask_gene()
    if gene is None:
        return
    conn = database()
    strand = {1: "+", -1: "-"}.get(gene["strand"], "?")
    print("\n" + Style.BRIGHT + tr("database.gene.header", symbol=gene["symbol"], name=gene["name"],
                                   organism=gene["scientific_name"]) + Style.RESET_ALL)
    if gene["start_pos"]:
        print("  " + tr("database.gene.location", chromosome=gene["chromosome"],
                        start=gene["start_pos"], end=gene["end_pos"], strand=strand,
                        assembly=gene["assembly"], gene_id=str(gene["ncbi_gene_id"])))

    pending = 0
    for seq in queries.sequences(conn, gene["id"]):
        if not seq["residues"]:
            pending += 1
            continue
        unit = tr("database.units.aa" if seq["mol_type"] == "protein" else "database.units.nt")
        line = tr("database.gene.sequence", accession=f"{seq['accession']}.{seq['version']}",
                  kind=tr(f"database.mol_types.{seq['mol_type']}"), length=seq["length"], unit=unit)
        if seq["cds_start"]:
            line += " · " + tr("database.gene.cds", start=seq["cds_start"], end=seq["cds_end"])
        print("\n  " + Fore.CYAN + line + Style.RESET_ALL)
        start = (seq["cds_start"] or 1) - 1
        print("    " + seq["residues"][start:start + LINE_WIDTH] + "…")
    if pending:
        print("\n  " + tr("database.gene.pending", count=pending))  # nº al final: sin plurales

    counts = queries.significance_counts(conn, gene["id"])
    if not counts:
        print("\n" + tr("database.gene.no_variants"))
        return
    print("\n" + tr("database.gene.variants", count=sum(n for _, n in counts)))
    for value, n in counts:
        print(f"  {tr.number(n):>7}  {significance(value)}")

    rows = queries.diseases(conn, gene["id"], DISEASES_SHOWN)
    if rows:
        print("\n" + tr("database.gene.diseases"))
        for row in rows:
            parts = [tr.count(row[effect], f"variant_{effect}")
                     for effect in ("causal", "protective", "risk_factor") if row[effect]]
            print(f"  {row['name_en']}: " + ", ".join(parts))
        print(Style.DIM + tr("database.gene.diseases_note") + Style.RESET_ALL)


def choose_variant(gene):
    """Busca variantes del gen por texto y deja elegir una. Devuelve la fila o None."""
    conn = database()
    text = input(tr("prompts.variant")).strip()
    found = queries.find_variants(conn, gene["id"], text, causal_only=not text, limit=VARIANTS_SHOWN)
    if not found:
        warning(tr("database.variants.none"))
        return None
    if len(found) == 1:
        return found[0]
    print("\n" + tr("database.variants.found"))
    for number, v in enumerate(found, start=1):
        name = f" «{v['common_name']}»" if v["common_name"] else ""
        print(f"  {number:>2}] {v['hgvs_c']} {v['hgvs_p'] or ''}{name} · {significance(v['clinical_significance'])}")
    if len(found) == VARIANTS_SHOWN:
        print(Style.DIM + tr("database.variants.more", count=VARIANTS_SHOWN) + Style.RESET_ALL)
    return found[ask_int(tr("prompts.choose_variant"), minimum=1, maximum=len(found)) - 1]


def db_compare():
    """La secuencia sana al lado de la variante, codón a codón, y lo que hace en cada enfermedad."""
    gene = ask_gene()
    if gene is None:
        return
    variant = choose_variant(gene)
    if variant is None:
        return
    conn = database()
    accession = f"{variant['accession']}.{variant['version']}"
    print("\n" + Style.BRIGHT + tr("database.compare.title", gene=gene["symbol"],
                                   hgvs_c=variant["hgvs_c"], hgvs_p=variant["hgvs_p"] or "",
                                   accession=accession) + Style.RESET_ALL)
    print("  " + tr("database.compare.global", significance=significance(variant["clinical_significance"])))
    try:
        NCBIProvider().ensure_residues(conn, variant["sequence_id"])
        effect, residues, cds_start = locate(conn, variant["id"])
    except DNAinventError as exc:
        warning(tr.error(exc))
    else:
        show_side_by_side(effect, side_by_side(effect, residues, cds_start))
    show_variant_diseases(variant)


def show_side_by_side(effect, view):
    print("\n  " + tr("database.compare.position", position=effect.position,
                      codon=effect.codon_number, ref=effect.ref_codon, alt=effect.alt_codon))
    labels = [tr("database.compare.codon"), tr("database.compare.healthy"), tr("database.compare.variant")]
    width = max(len(label) for label in labels) + 2
    numbers = "".join(f"{view.first_codon + i:<5}" for i in range(len(view.healthy)))
    print(f"\n  {labels[0]:<{width}}" + Style.DIM + numbers + Style.RESET_ALL)
    for label, codons, protein in ((labels[1], view.healthy, view.healthy_protein),
                                   (labels[2], view.variant, view.variant_protein)):
        line, aas = [], []
        for i, codon in enumerate(codons):
            if i == view.changed_codon:
                b = view.changed_base
                color = Fore.GREEN if codons is view.healthy else Fore.RED
                codon = (Fore.CYAN + codon[:b] + color + Style.BRIGHT + codon[b] + Style.RESET_ALL
                         + Fore.CYAN + codon[b + 1:] + Style.RESET_ALL)
                aa = Style.BRIGHT + three_letter(protein[i]) + Style.RESET_ALL
            else:
                aa = Style.DIM + three_letter(protein[i]) + Style.RESET_ALL
            line.append(codon + "  ")
            aas.append(aa + "  ")
        print(f"  {label:<{width}}" + "".join(line))
        print(f"  {'':<{width}}" + "".join(aas))
    print("\n  " + tr("database.compare.effect", kind=tr(f"database.variant_kinds.{effect.kind}")))
    if effect.kind == "missense" and effect.codon_number > 1:
        old = f"{three_letter(effect.ref_aa)}{effect.codon_number - 1}{three_letter(effect.alt_aa)}"
        print(Style.DIM + "  " + tr("database.compare.old_numbering", old=old) + Style.RESET_ALL)


def show_variant_diseases(variant):
    rows = queries.variant_diseases(database(), variant["id"])
    if not rows:
        print("\n" + tr("database.compare.no_diseases"))
        return
    print("\n" + tr("database.compare.diseases"))
    for row in rows:
        effect = tr(f"database.effects.{row['effect'] or 'none'}")
        print(f"  - {row['name_en']}: {effect} ({significance(row['clinical_significance'])})")


def db_search():
    text = input(tr("prompts.search_text")).strip()
    if not text:
        return
    conn = database()
    rows = queries.search(conn, text)
    if not rows:
        warning(tr("database.search.none"))
        return
    for row in rows:
        if row["kind"] == "gene":
            print("  " + Fore.CYAN + tr("database.search.gene") + Style.RESET_ALL + " " + row["label"])
            continue
        print("  " + Fore.MAGENTA + tr("database.search.disease") + Style.RESET_ALL + " " + row["label"])
        genes = [f"{g['symbol']} ({tr.count(g['causal'], 'variant_causal')})"
                 for g in queries.genes_for_disease(conn, row["ref_id"])]
        if genes:
            print("      " + tr("database.search.genes", list=", ".join(genes)))


def db_list():
    rows = queries.genes(database())
    if not rows:
        warning(tr("database.empty"))
    for row in rows:
        print(f"  {Fore.CYAN}{row['symbol']:<10}{Style.RESET_ALL} {row['name']} · "
              f"{row['scientific_name']} · {tr.count(row['variant_count'], 'variant')}")
    print(Style.DIM + tr("database.path", path=default_path()) + Style.RESET_ALL)


def database_menu():
    menu("app.database_title", {
        "1": ("menu.database.import", db_import),
        "2": ("menu.database.gene", db_gene),
        "3": ("menu.database.compare", db_compare),
        "4": ("menu.database.search", db_search),
        "5": ("menu.database.list", db_list),
    }, back_key="menu.back")


def option_language():
    choose_language()
    print(tr("language.changed"))


def choose_language():
    """Muestra los idiomas, cambia al elegido y lo guarda para la próxima vez."""
    languages = available_languages()
    print("\n" + " / ".join(Translator(code)("language.choose") for code in languages))
    codes = list(languages)
    for number, code in enumerate(codes, start=1):
        print(f"    {number}] {languages[code]}")
    choice = ask_int(tr("common.choose_option"), minimum=1, maximum=len(codes))
    set_language(codes[choice - 1])
    settings.save({"language": tr.language})


def set_language(language: str):
    global tr
    tr = Translator(language)


# ---------- Submenú de herramientas ----------

def tool_find():
    dna = ask_dna()
    target = ask_dna("prompts.search")
    positions = tools.find_positions(dna, target)
    print(highlight(dna, positions, len(target)))
    print(tr("find.found", count=len(positions), positions=[p + 1 for p in positions]))


def tool_similarity():
    first = ask_dna("prompts.first_dna")
    second = ask_dna("prompts.second_dna")
    cells = len(first) * len(second)
    if cells > ALIGNMENT_WARNING and not confirm("align.confirm", cells=cells):
        return

    result = alignment.align(first, second)
    show_alignment(result)
    print_labeled([
        ("align.identity", f"{tr.decimal(result.identity, 2)} % ({result.matches}/{result.length})"),
        ("align.mismatches", str(result.mismatches)),
        ("align.gaps", str(result.gaps)),
        ("align.score", str(result.score)),
    ])


def show_alignment(result: alignment.Alignment):
    print()
    for start in range(0, result.length, LINE_WIDTH):
        end = start + LINE_WIDTH
        top, mid, bottom = result.first[start:end], result.markup[start:end], result.second[start:end]
        print("  " + "".join(color_aligned(a, b) for a, b in zip(top, bottom)) + Style.RESET_ALL)
        print("  " + mid)
        print("  " + "".join(color_aligned(b, a) for a, b in zip(top, bottom)) + Style.RESET_ALL)
        print()


def color_aligned(base: str, other: str) -> str:
    """Verde si coincide, amarillo si es una sustitución, rojo si hay hueco."""
    if alignment.GAP in (base, other):
        return Fore.RED + base
    return (Fore.GREEN if base == other else Fore.YELLOW) + base


def tool_permutations():
    dna = ask_dna()
    r = ask_int(tr("prompts.permutation_size"), minimum=1)
    if r > len(dna):
        print(tr("permutations.too_big", size=r, length=len(dna)))
        return
    total = tools.count_permutations(dna, r)
    positional = tools.count_positional_permutations(dna, r)
    print_labeled([("permutations.positional", tr.number(positional)),
                   ("permutations.distinct", tr.number(total))])
    if positional > total:
        print(tr("permutations.avoided", count=positional - total))
    if total > PERMUTATION_WARNING and not confirm("permutations.confirm"):
        return
    for p in tools.iter_permutations(dna, r):
        print(p)


def show_being(title_key: str, bit_string: str):
    print("\n" + tr(title_key))
    print_labeled([("bits.label_bits", color_bits(bit_string)),
                   ("bits.label_dna", color_dna(sequences.from_binary(bit_string)))], indent="  ")


def tool_bits():
    n = ask_int(tr("prompts.being_size"), minimum=1)
    original, _, elapsed = bits.generate_being(n)
    current = original  # estado actual: cada operación parte de aquí y lo actualiza
    show_being("bits.generated", current)

    def apply(title_key: str, transform) -> str:
        """Aplica transform al estado actual, lo guarda y devuelve el estado anterior."""
        nonlocal current
        previous, current = current, transform(current)
        show_being(title_key, current)
        return previous

    def reverse_complement_step():
        previous = apply("bits.both", lambda b: bits.flip(bits.reverse_blocks(b)))
        check = sequences.from_binary(current) == sequences.reverse_complement(
            sequences.from_binary(previous))
        print("  " + tr("bits.check", result=check))

    def restore():
        nonlocal current
        current = original
        show_being("bits.restored", current)

    def show_time():
        print("\n" + tr("bits.time", bits=len(original),
                        ms=tr.decimal(elapsed * 1000, 3),
                        us=tr.decimal(elapsed / len(original) * 1e6, 3)))

    menu("app.bits_title", {
        "1": ("menu.bits.time", show_time),
        "2": ("menu.bits.flip", lambda: apply("bits.flip", bits.flip)),
        "3": ("menu.bits.reverse", lambda: apply("bits.reverse", bits.reverse_blocks)),
        "4": ("menu.bits.both", reverse_complement_step),
        "5": ("menu.bits.restore", restore),
    }, back_key="menu.back")


def tools_menu():
    menu("app.tools_title", {
        "1": ("menu.tools.find", tool_find),
        "2": ("menu.tools.align", tool_similarity),
        "3": ("menu.tools.permutations", tool_permutations),
        "4": ("menu.tools.bits", tool_bits),
    }, back_key="menu.back")


# ---------- Motor de menús ----------

def menu(title_key: str, items: dict[str, tuple[str, callable]], back_key: str) -> None:
    """Muestra un menú hasta que el usuario elige 0.

    items: tecla -> (clave del texto, función). Los textos se traducen cada vez
    que se dibuja el menú, así un cambio de idioma se ve al momento.
    """
    while True:
        print(f"\n        {tr(title_key)}")
        for key, (label_key, _) in items.items():
            print(f"    {key}] {tr(label_key)}")
        print(f"    0] {tr(back_key)}")
        choice = input(tr("common.choose_option")).strip().upper()
        if choice == "0":
            return
        if choice not in items:
            print(tr("common.invalid_option"))
            continue
        items[choice][1]()


def run(language: str | None = None) -> None:
    """language: el de --lang. Si no se indica, el guardado; y la primera vez, se pregunta."""
    init()
    try:
        saved = settings.language()
        if language:
            set_language(language)
        elif saved in available_languages():
            set_language(saved)
        else:
            set_language(DEFAULT_LANGUAGE)
            choose_language()
        menu("app.title", {
            "1": ("menu.main.complement", option_complement),
            "2": ("menu.main.transcribe", option_transcribe),
            "3": ("menu.main.translate", option_translate),
            "4": ("menu.main.crispr", option_crispr),
            "5": ("menu.main.to_binary", option_binary),
            "6": ("menu.main.from_binary", option_from_binary),
            "7": ("menu.main.microscope", option_microscope),
            "8": ("menu.main.protein", option_protein),
            "9": ("menu.main.database", database_menu),
            "*": ("menu.main.tools", tools_menu),
            "L": ("menu.main.language", option_language),
        }, back_key="menu.main.exit")
    except (KeyboardInterrupt, EOFError):
        print()
    print(tr("app.goodbye"))
    sys.exit(0)
