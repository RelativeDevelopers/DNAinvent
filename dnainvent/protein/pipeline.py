"""Del ADN del usuario a su proteína: el recorrido completo, sin interfaz.

    report = analyze("ATGGTGCATCTG...", blast=BlastClient(email), afdb=AlphaFoldDB())

    1. ORF        el marco de lectura más largo de las 6 pautas, traducido     (orf.py)
    2. BLAST      parecidos en Swiss-Prot; si no hay ninguno fiable, en todo   (clients/blast.py)
                  UniProtKB
    3. Caso       conocida / variante / pariente lejano / inventada            (classify.py)
    4. Estructura la de AlphaFold DB de la proteína real más parecida          (clients/afdb.py)
    5. Predicción la de la proteína del usuario, si no es idéntica a una real  (predictors/)
                  (variante, pariente lejano o inventada)
    6. Comparación mutaciones y diferencia de forma con la natural            (compare.py)

Como todo el paquete, no imprime nada: devuelve un ProteinReport y, mientras trabaja, avisa
de cada paso con progress(clave, datos). Las claves son de i18n (protein.stage.*), así cada
interfaz lo cuenta en su idioma. Los avisos para el usuario también son claves
(protein.warnings.*) con sus datos.

Que falle la estructura o la predicción no estropea el resto (se avisa y el informe sigue).
Que falle BLAST sí: sin él no se sabe qué proteína es, así que el error sube a la interfaz.
"""

from dataclasses import dataclass, field
from typing import Callable

from ..errors import DNAinventError, ServiceError
from . import confidence
from .classify import Category, Classification, Thresholds, classify
from .clients.afdb import AlphaFoldDB, AlphaFoldEntry
from .clients.blast import SWISSPROT, UNIPROTKB, BlastClient, BlastHit
from .compare import Comparison, compare
from .orf import ORF, best_orf
from .pdb import Structure
from .predictors.base import Predictor

# Menos de esto no se puede identificar: un parecido perfecto de 15 aminoácidos todavía
# podría ser casualidad (su E-value no baja de ~1e-3). Tampoco se predice su forma: un
# péptido tan corto no tiene un plegamiento propio (y ESMFold suele rechazarlo).
MIN_SEARCH_LENGTH = 20
UNLIKELY_FOLD = 50.0      # pLDDT medio por debajo del cual "probablemente no se plegaría"

Progress = Callable[[str, dict], None]


@dataclass
class ProteinReport:
    dna_length: int
    orf: ORF
    classification: Classification
    database: str | None = None              # dónde se encontró (o la última donde se buscó)
    hits: list[BlastHit] = field(default_factory=list)
    best: BlastHit | None = None              # el parecido que decide el caso
    entry: AlphaFoldEntry | None = None
    structure: Structure | None = None        # la natural (AlphaFold DB)
    prediction: Structure | None = None       # la del usuario (predicha)
    predictor: str | None = None              # quién la predijo ("ESMFold")
    comparison: Comparison | None = None      # tuya frente a la natural
    warnings: list[tuple[str, dict]] = field(default_factory=list)   # (clave i18n, datos)

    @property
    def protein(self) -> str:
        return self.orf.protein

    @property
    def category(self) -> Category:
        return self.classification.category


def analyze(dna: str, blast: BlastClient, afdb: AlphaFoldDB | None = None,
            predictor: Predictor | None = None, progress: Progress | None = None,
            thresholds: Thresholds = Thresholds(),
            search_everything: bool = True) -> ProteinReport:
    """Recorrido completo. search_everything: si Swiss-Prot no da nada fiable, buscar
    también en UniProtKB completa (más lenta). Sin predictor no se predice nada."""
    say = progress or (lambda key, data: None)

    say("protein.stage.orf", {})
    orf = best_orf(dna)                       # InvalidSequenceError("no_orf") si no hay ATG
    report = ProteinReport(len(dna), orf, Classification(Category.INVENTED, None))
    if not orf.is_complete:
        report.warnings.append(("protein.warnings.incomplete_orf", {}))

    if orf.length < MIN_SEARCH_LENGTH:
        report.warnings.append(("protein.warnings.too_short",
                                {"length": orf.length, "minimum": MIN_SEARCH_LENGTH}))
        _category_warnings(report)
        return report

    databases = (SWISSPROT, UNIPROTKB) if search_everything else (SWISSPROT,)
    for database in databases:
        say("protein.stage.search" if database == SWISSPROT else "protein.stage.search_everything",
            {"database": database})
        result = blast.search(orf.protein, database,
                              progress=lambda status: say("protein.stage.waiting",
                                                          {"database": database,
                                                           "status": status}))
        report.database, report.hits = database, result.hits
        report.classification = classify([hit.to_hit() for hit in result.hits], thresholds)
        if report.classification.is_known:
            break

    if report.classification.is_known:
        accession = report.classification.hit.accession
        report.best = next(h for h in report.hits if h.accession == accession)
        if report.best.subject_coverage < 100 - 1e-9:
            report.warnings.append(("protein.warnings.fragment", {
                "start": report.best.hit_from, "end": report.best.hit_to,
                "length": report.best.length}))
        if afdb is not None:
            say("protein.stage.structure", {"accession": accession})
            try:
                report.entry, report.structure = afdb.structure(accession)
            except ServiceError as error:
                key = ("protein.warnings.no_structure" if error.code == "protein_not_found"
                       else "protein.warnings.structure_unavailable")
                report.warnings.append((key, {"accession": accession}))

    if report.category is not Category.IDENTICAL:
        _predict(report, predictor, say)
    if report.best is not None and report.category is not Category.IDENTICAL:
        report.comparison = compare(report.protein, _natural_sequence(report), report.best,
                                    report.structure, report.prediction)
    _category_warnings(report)
    return report


def _predict(report: ProteinReport, predictor: Predictor | None, say: Progress) -> None:
    if predictor is None:
        return
    say("protein.stage.predict", {"predictor": predictor.name, "length": report.orf.length})
    try:
        report.prediction = predictor.predict(report.protein)
        report.predictor = predictor.name
    except DNAinventError as error:
        if error.code == "prediction_too_long":
            report.warnings.append(("protein.warnings.prediction_too_long", dict(error.params)))
        elif error.code == "prediction_failed":
            report.warnings.append(("protein.warnings.prediction_failed",
                                    {"predictor": predictor.name}))
        else:
            report.warnings.append(("protein.warnings.prediction_unavailable",
                                    {"predictor": predictor.name}))


def _natural_sequence(report: ProteinReport) -> str:
    """Secuencia completa de la natural. Si no se pudo descargar su estructura, se rehace
    con lo que dice BLAST (las posiciones que no salen en el alineamiento quedan como X)."""
    if report.entry is not None and report.entry.sequence:
        return report.entry.sequence
    if report.structure is not None:
        return report.structure.sequence()
    hit = report.best
    letters = ["X"] * hit.length
    position = hit.hit_from - 1
    for letter in hit.hit_aligned:
        if letter != "-":
            letters[position] = letter
            position += 1
    return "".join(letters)


def _category_warnings(report: ProteinReport) -> None:
    if report.category is Category.VERY_SIMILAR:
        report.warnings.append(("protein.warnings.mutation_effect", {}))
    if report.category is Category.INVENTED and report.prediction is not None:
        plddt = confidence.mean(report.prediction)
        if plddt < UNLIKELY_FOLD:
            report.warnings.append(("protein.warnings.unlikely_fold",
                                    {"plddt": round(plddt, 1)}))
