"""Los genes de serie: la base de datos que DNAinvent trae ya hecha (consola y web).

    builder = CatalogBuilder(open_database("catalog.db"))
    report = builder.run()                         # los 14 genes; se puede cortar y seguir

Se construye en el ordenador de quien publica (tools/build_catalog.py), porque hace miles de
consultas al NCBI y a ChEMBL. La web solo la LEE: así el servidor no llama a ninguna API.

Pasos de cada gen, en orden:
    gene           NCBI Gene: el gen, su ARNm de referencia y su proteína
    variants       ClinVar: todas sus variantes, con la clasificación por enfermedad y la
                   frecuencia en la población
    complete       vuelve a buscar en ClinVar y recupera las variantes del gen que falten
                   (resúmenes que no llegaron); también arregla catálogos antiguos
    sequences      las secuencias de los ARNm sobre los que se describen las variantes (para
                   poder comparar sana frente a variante sin conexión)
    gene_drugs     ChEMBL: fármacos que actúan sobre su proteína
    disease_drugs  ChEMBL: fármacos para sus enfermedades principales
Y al final, para todo el catálogo:
    combined       quita los registros de ClinVar que son combinaciones de variantes y se
                   guardaron como una sola (catálogos hechos antes de filtrarlas)
    orphanet       herencia, edad de inicio y prevalencia de las enfermedades raras
    common_names   nombres conocidos de algunas variantes (HbS, ΔF508...)

Cada paso terminado se apunta en la tabla catalog_steps: si se corta (sin conexión, ChEMBL
caído...), al volver a ejecutarlo continúa donde lo dejó. Un paso que falla no se apunta y no
para los demás.
"""

import sqlite3
from dataclasses import dataclass, field
from typing import Callable

from .. import __version__
from ..errors import DNAinventError
from . import queries
from .chembl import ChEMBLProvider
from .clinvar import ClinVarProvider, is_combined
from .ncbi import HUMAN, NCBIProvider
from .orphanet import OrphanetProvider
from .records import now_iso

Progress = Callable[[str, dict], None]

CATALOG_GENES = ("HBB", "CFTR", "BRCA1", "BRCA2", "TP53", "HTT", "DMD", "F8", "PAH", "APOE",
                 "LDLR", "SMN1", "G6PD", "LCT")
GENE_STEPS = ("gene", "variants", "complete", "sequences", "gene_drugs", "disease_drugs")
DISEASES_WITH_DRUGS = 5       # enfermedades de cada gen para las que se buscan fármacos

# Nombres conocidos de variantes: (gen, HGVS c., HGVS p., nombre). Se buscan primero por c. y,
# si no está, por p. Los nombres son los habituales en la literatura; su explicación en cada
# idioma va aparte (i18n), no aquí.
#   HbS / HbC: hemoglobina de la anemia falciforme y la hemoglobina C (codón 6 en la numeración
#              antigua, sin la metionina inicial; 7 en HGVS).
#   185delAG, 5382insC (BRCA1) y 6174delT (BRCA2): las tres variantes fundadoras asquenazíes,
#              con su nombre antiguo (numeración de BIC).
#   G6PD Mediterranean: ClinVar la nombra sobre el transcrito antiguo NM_000402.4 (c.653C>T,
#              p.Ser218Phe); en el actual es c.563C>T (p.Ser188Phe). Es la misma variante
#              (rs5030868), así que se aceptan las dos formas.
COMMON_NAMES = (
    ("HBB", "c.20A>T", "p.Glu7Val", "HbS"),
    ("HBB", "c.19G>A", "p.Glu7Lys", "HbC"),
    ("CFTR", "c.1521_1523del", "p.Phe508del", "ΔF508 (F508del)"),
    ("CFTR", "c.1652G>A", "p.Gly551Asp", "G551D"),
    ("BRCA1", "c.68_69del", "p.Glu23fs", "185delAG"),
    ("BRCA1", "c.5266dup", "p.Gln1756fs", "5382insC"),
    ("BRCA2", "c.5946del", "p.Ser1982fs", "6174delT"),
    ("TP53", "c.524G>A", "p.Arg175His", "R175H"),
    ("TP53", "c.818G>A", "p.Arg273His", "R273H"),
    ("PAH", "c.1222C>T", "p.Arg408Trp", "R408W"),
    ("APOE", "c.388T>C", "p.Cys130Arg", "APOE ε4"),
    ("APOE", "c.526C>T", "p.Arg176Cys", "APOE ε2"),
    ("G6PD", "c.653C>T", "p.Ser218Phe", "G6PD Mediterranean"),
    ("G6PD", "c.563C>T", "p.Ser188Phe", "G6PD Mediterranean"),
)

_TABLES = """
CREATE TABLE IF NOT EXISTS catalog_steps (
    symbol      TEXT NOT NULL,          -- gen, o '*' para los pasos de todo el catálogo
    step        TEXT NOT NULL,
    finished_at TEXT NOT NULL,
    detail      TEXT,
    PRIMARY KEY (symbol, step)
);
CREATE TABLE IF NOT EXISTS catalog_info (
    key   TEXT PRIMARY KEY,
    value TEXT
);
"""


def _no_progress(key: str, data: dict) -> None:
    pass


def _belongs_to(description: str | None, symbol: str) -> bool:
    """¿La definición de RefSeq es de ese gen? Llevan el símbolo entre paréntesis:
    «Homo sapiens huntingtin (HTT), transcript variant 2, mRNA»."""
    return description is None or f"({symbol.upper()})" in description.upper()


class StepIncomplete(Exception):
    """Un paso que ha ido a medias: no se apunta como hecho, para repetirlo la próxima vez.
    El texto es para quien construye el catálogo (herramienta de desarrollo, en castellano)."""


@dataclass
class StepResult:
    symbol: str
    step: str
    status: str               # 'done', 'skipped' (ya estaba), 'failed'
    detail: str = ""


@dataclass
class CatalogReport:
    steps: list[StepResult] = field(default_factory=list)
    missing_names: list[tuple[str, str]] = field(default_factory=list)   # (gen, nombre) no encontrados
    removed_genes: list[str] = field(default_factory=list)              # genes borrados por limpieza
    repaired: list[str] = field(default_factory=list)                   # ARNm de referencia corregidos

    @property
    def failed(self) -> list[StepResult]:
        return [s for s in self.steps if s.status == "failed"]


def catalog_genes(conn: sqlite3.Connection) -> list[str]:
    """Genes que ya tienen todos sus pasos hechos en esta base de datos."""
    try:
        rows = conn.execute("SELECT symbol, COUNT(*) AS n FROM catalog_steps WHERE symbol != '*' "
                            "GROUP BY symbol").fetchall()
    except sqlite3.OperationalError:          # una base de datos normal, sin catálogo
        return []
    return [r["symbol"] for r in rows if r["n"] == len(GENE_STEPS)]


class CatalogBuilder:
    def __init__(self, conn: sqlite3.Connection, ncbi: NCBIProvider | None = None,
                 clinvar: ClinVarProvider | None = None, chembl: ChEMBLProvider | None = None,
                 orphanet: OrphanetProvider | None = None, progress: Progress = _no_progress,
                 max_variants: int | None = None, diseases_with_drugs: int = DISEASES_WITH_DRUGS,
                 announce: Callable[[StepResult], None] | None = None,
                 on_gene: Callable[[str], None] | None = None):
        self.conn = conn
        self.ncbi = ncbi or NCBIProvider()
        self.clinvar = clinvar or ClinVarProvider(http=self.ncbi.http, email=self.ncbi.email,
                                                  api_key=self.ncbi.api_key, sleep=self.ncbi.sleep)
        self.chembl = chembl or ChEMBLProvider()
        self.orphanet = orphanet or OrphanetProvider()
        self.progress = progress
        self.max_variants = max_variants
        self.diseases_with_drugs = diseases_with_drugs
        self.announce = announce or (lambda result: None)
        self.on_gene = on_gene or (lambda symbol: None)     # antes de empezar cada gen ('*': el final)
        self.repaired: list[str] = []      # genes con un ARNm de referencia equivocado, ya quitado
        conn.executescript(_TABLES)

    # ------------------------------------------------------------------ pasos
    def done(self, symbol: str, step: str) -> bool:
        return self.conn.execute("SELECT 1 FROM catalog_steps WHERE symbol = ? AND step = ?",
                                 (symbol, step)).fetchone() is not None

    def _mark(self, symbol: str, step: str, detail: str = "") -> None:
        with self.conn:
            self.conn.execute("INSERT OR REPLACE INTO catalog_steps VALUES (?, ?, ?, ?)",
                              (symbol, step, now_iso(), detail))

    def _gene_id(self, symbol: str) -> int | None:
        row = self.conn.execute("SELECT id FROM genes WHERE symbol = ? COLLATE NOCASE AND taxon_id = ?",
                                (symbol, HUMAN)).fetchone()
        return row["id"] if row else None

    # ------------------------------------------------------------------ limpieza
    def forget_gene(self, gene_id: int) -> None:
        """Borra un gen y todo lo suyo (variantes, secuencias, dianas...) y las enfermedades que
        se queden sin nada. Para genes que no deberían estar en el catálogo."""
        c = self.conn
        with c:
            seqs = [r[0] for r in c.execute("SELECT id FROM sequences WHERE gene_id = ?", (gene_id,))]
            c.execute("DELETE FROM variant_diseases WHERE variant_id IN "
                      "(SELECT id FROM variants WHERE gene_id = ?)", (gene_id,))
            c.execute("DELETE FROM variants WHERE gene_id = ?", (gene_id,))
            c.execute("DELETE FROM gene_diseases WHERE gene_id = ?", (gene_id,))
            c.execute("DELETE FROM drug_targets WHERE gene_id = ?", (gene_id,))
            for seq in seqs:
                c.execute("DELETE FROM sequence_links WHERE from_sequence_id = ? OR to_sequence_id = ?",
                          (seq, seq))
                c.execute("DELETE FROM sequence_tags WHERE sequence_id = ?", (seq,))
                c.execute("UPDATE drugs SET sequence_id = NULL WHERE sequence_id = ?", (seq,))
            c.execute("DELETE FROM sequences WHERE gene_id = ?", (gene_id,))
            c.execute("DELETE FROM search_index WHERE kind = 'gene' AND ref_id = ?", (gene_id,))
            c.execute("DELETE FROM genes WHERE id = ?", (gene_id,))
            orphans = [r[0] for r in c.execute(
                """SELECT id FROM diseases WHERE id NOT IN (SELECT disease_id FROM variant_diseases)
                   AND id NOT IN (SELECT disease_id FROM gene_diseases)
                   AND id NOT IN (SELECT disease_id FROM drug_diseases)""")]
            for disease_id in orphans:
                c.execute("DELETE FROM search_index WHERE kind = 'disease' AND ref_id = ?", (disease_id,))
                c.execute("DELETE FROM diseases WHERE id = ?", (disease_id,))

    def clean(self) -> list[str]:
        """Deja la base de datos coherente antes de seguir:
        - borra los genes que no son del catálogo (p. ej. SLC6A4, que entró por error al buscar
          HTT, cuando la búsqueda aceptaba alias);
        - olvida los pasos hechos de un gen que ya no está.
        Devuelve los símbolos borrados."""
        wanted = {g.upper() for g in CATALOG_GENES}
        removed = []
        for row in self.conn.execute("SELECT id, symbol FROM genes").fetchall():
            if (row["symbol"] or "").upper() not in wanted:
                self.forget_gene(row["id"])
                removed.append(row["symbol"])
        with self.conn:
            for (symbol,) in self.conn.execute(
                    "SELECT DISTINCT symbol FROM catalog_steps WHERE symbol != '*'").fetchall():
                if self._gene_id(symbol) is None:
                    self.conn.execute("DELETE FROM catalog_steps WHERE symbol = ?", (symbol,))
        # Un ARNm de referencia que no es del gen (HTT tenía el de SLC6A4): se quita y se
        # repiten el gen y la búsqueda de fármacos por su proteína. Las variantes no cambian.
        for row in self.conn.execute("SELECT id, symbol FROM genes").fetchall():
            mrna = self._reference_mrna(row["id"])
            if mrna is not None and not _belongs_to(mrna["description"], row["symbol"]):
                self._drop_reference(row["id"])
                with self.conn:
                    self.conn.execute("DELETE FROM catalog_steps WHERE symbol = ? AND step IN "
                                      "('gene', 'gene_drugs')", (row["symbol"],))
                self.repaired.append(f"{row['symbol']} ({mrna['accession']}.{mrna['version']})")
        return removed

    def _step(self, report: CatalogReport, symbol: str, step: str, action, redo: bool) -> bool:
        """Hace un paso (si no estaba hecho) y lo apunta. False si ha fallado."""
        if not redo and self.done(symbol, step):
            result = StepResult(symbol, step, "skipped")
        else:
            try:
                detail = action() or ""
            except StepIncomplete as exc:
                result = StepResult(symbol, step, "failed", str(exc))
            except DNAinventError as exc:
                result = StepResult(symbol, step, "failed", str(exc))
            except Exception as exc:          # un fallo del programa no debe parar todo el catálogo
                result = StepResult(symbol, step, "failed", f"error interno: {type(exc).__name__}: {exc}")
            else:
                self._mark(symbol, step, detail)
                result = StepResult(symbol, step, "done", detail)
        report.steps.append(result)
        self.announce(result)
        return result.status != "failed"

    def _reference_mrna(self, gene_id: int):
        return self.conn.execute("SELECT id, accession, version, description FROM sequences "
                                 "WHERE gene_id = ? AND is_reference = 1 AND mol_type = 'mrna'",
                                 (gene_id,)).fetchone()

    def _import_gene(self, symbol: str) -> str:
        gene_id = self._gene_id(symbol)
        if gene_id is None or self._reference_mrna(gene_id) is None:
            gene_id = self.ncbi.import_gene(self.conn, symbol, HUMAN, self.progress)
        row = self.conn.execute("SELECT symbol, name FROM genes WHERE id = ?", (gene_id,)).fetchone()
        if (row["symbol"] or "").upper() != symbol.upper():
            found = row["symbol"]
            self.forget_gene(gene_id)
            raise StepIncomplete(f"el NCBI ha devuelto {found} en lugar de {symbol}")
        mrna = self._reference_mrna(gene_id)
        if mrna is not None and not _belongs_to(mrna["description"], symbol):
            self._drop_reference(gene_id)
            raise StepIncomplete(f"el ARNm de referencia {mrna['accession']}.{mrna['version']} "
                                 f"no es de {symbol}: «{mrna['description']}»")
        detail = row["name"] or ""
        return f"{detail} · {mrna['accession']}.{mrna['version']}" if mrna else detail

    def _drop_reference(self, gene_id: int) -> None:
        """Quita un ARNm de referencia equivocado y su proteína (las variantes no lo usan: se
        describen sobre los ARNm del gen de verdad). Si alguna lo usara, solo se desmarca."""
        c = self.conn
        with c:
            for seq in c.execute("SELECT id FROM sequences WHERE gene_id = ? AND is_reference = 1",
                                 (gene_id,)).fetchall():
                used = c.execute("SELECT 1 FROM variants WHERE sequence_id = ? LIMIT 1",
                                 (seq["id"],)).fetchone()
                if used:
                    c.execute("UPDATE sequences SET is_reference = 0 WHERE id = ?", (seq["id"],))
                    continue
                c.execute("DELETE FROM sequence_links WHERE from_sequence_id = ? OR to_sequence_id = ?",
                          (seq["id"], seq["id"]))
                c.execute("DELETE FROM sequence_tags WHERE sequence_id = ?", (seq["id"],))
                c.execute("UPDATE drugs SET sequence_id = NULL WHERE sequence_id = ?", (seq["id"],))
                c.execute("DELETE FROM sequences WHERE id = ?", (seq["id"],))

    @staticmethod
    def _describe(stats) -> str:
        """«838 guardadas de 1180 (300 sin ARNm, 40 combinaciones...); 343 con frecuencia»"""
        reasons = [(stats.not_on_mrna, "sin ARNm"), (stats.combined, "combinaciones"),
                   (stats.other_gene, "de otro gen"), (stats.summaries_missing, "sin respuesta")]
        why = ", ".join(f"{n} {label}" for n, label in reasons if n)
        return (f"{stats.saved} guardadas de {stats.found}" + (f" ({why})" if why else "")
                + f"; {stats.with_frequency} con frecuencia")

    def _check_import(self, stats) -> None:
        if stats.summaries_missing:
            raise StepIncomplete(f"{self._describe(stats)}; faltan {stats.summaries_missing} "
                                 "resúmenes de ClinVar")

    def _variants(self, symbol: str) -> str:
        stats = self.clinvar.import_variants(self.conn, symbol, max_variants=self.max_variants,
                                             progress=self.progress)
        detail = self._describe(stats)
        self._check_import(stats)
        if stats.condition_batches_failed:
            # Sin la clasificación por enfermedad el gen queda cojo: que se repita
            raise StepIncomplete(f"{stats.condition_batches_failed} lotes de RCV fallidos")
        return detail

    def _complete(self, symbol: str) -> str:
        """Pide otra vez a ClinVar la lista de variantes del gen e importa las que no están
        guardadas. Casi todas resultan ser de las que no se guardan nunca (sin ARNm,
        combinaciones); las demás son las que se perdieron."""
        gene_id = self._gene_id(symbol)
        ids = self.clinvar.search_variant_ids(symbol, self.max_variants)
        have = {str(r[0]) for r in self.conn.execute(
            "SELECT clinvar_id FROM variants WHERE gene_id = ? AND clinvar_id IS NOT NULL", (gene_id,))}
        missing = [uid for uid in ids if uid not in have]
        if not missing:
            return "no falta ninguna"
        stats = self.clinvar.import_variants(self.conn, symbol, progress=self.progress,
                                             only_ids=missing)
        self._check_import(stats)
        if stats.condition_batches_failed:
            raise StepIncomplete(f"{stats.condition_batches_failed} lotes de RCV fallidos")
        return f"{stats.saved} recuperadas de {len(missing)} revisadas"

    def _sequences(self, symbol: str) -> str:
        gene_id = self._gene_id(symbol)
        rows = self.conn.execute(
            """SELECT DISTINCT s.id, s.accession, s.version FROM variants v
               JOIN sequences s ON s.id = v.sequence_id
               WHERE v.gene_id = ? AND s.residues IS NULL""", (gene_id,)).fetchall()
        failed = []
        for row in rows:
            self.progress("database.stage.sequence", {"accession": f"{row['accession']}.{row['version']}"})
            try:
                self.ncbi.ensure_residues(self.conn, row["id"])
            except DNAinventError:
                failed.append(f"{row['accession']}.{row['version']}")
        if failed:
            raise StepIncomplete("sin secuencia: " + ", ".join(failed))
        return f"{len(rows)} descargadas"

    def _gene_drugs(self, symbol: str) -> str:
        stats = self.chembl.import_gene_drugs(self.conn, self._gene_id(symbol), self.progress)
        return f"proteína {stats.accession or '-'}; {stats.drugs} fármacos"

    def _disease_drugs(self, symbol: str) -> str:
        diseases = [d for d in queries.diseases(self.conn, self._gene_id(symbol),
                                                limit=self.diseases_with_drugs) if d["causal"]]
        failed = 0
        for disease in diseases:
            try:
                self.chembl.import_disease_drugs(self.conn, disease["id"], self.progress)
            except DNAinventError:
                failed += 1
        if failed:
            raise StepIncomplete(f"{failed} de {len(diseases)} enfermedades sin respuesta de ChEMBL")
        return f"{len(diseases)} enfermedades"

    def _remove_combined(self) -> str:
        """Los catálogos hechos antes de filtrar las combinaciones (clinvar.is_combined) las
        tienen guardadas como variantes sueltas. Todas quedaron con el tipo 'other': se vuelve
        a pedir su resumen (200 por consulta) y se borran las que lo son."""
        rows = self.conn.execute("SELECT id, clinvar_id FROM variants WHERE variant_type = 'other' "
                                 "AND clinvar_id IS NOT NULL").fetchall()
        by_clinvar = {str(r["clinvar_id"]): r["id"] for r in rows}
        docs = self.clinvar.fetch_summaries(list(by_clinvar), self.progress)
        combined = [by_clinvar[uid] for uid, doc in docs if uid in by_clinvar and is_combined(doc)]
        with self.conn:
            for variant_id in combined:
                self.conn.execute("DELETE FROM variant_diseases WHERE variant_id = ?", (variant_id,))
                self.conn.execute("DELETE FROM variants WHERE id = ?", (variant_id,))
        if len(docs) < len(by_clinvar):
            raise StepIncomplete(f"ClinVar solo ha respondido {len(docs)} de {len(by_clinvar)}; "
                                 f"borradas {len(combined)}")
        return f"{len(combined)} combinaciones borradas de {len(by_clinvar)} registros revisados"

    def _orphanet(self) -> str:
        stats = self.orphanet.apply(self.conn, self.progress)
        with self.conn:
            for kind, date in stats.file_dates.items():
                self.conn.execute("INSERT OR REPLACE INTO catalog_info VALUES (?, ?)",
                                  (f"orphadata_{kind}_date", date))
        return (f"{stats.by_code + stats.by_name} de {stats.diseases} enfermedades en Orphanet; "
                f"{stats.inheritance} con herencia, {stats.prevalence} con prevalencia")

    def apply_common_names(self, report: CatalogReport | None = None) -> tuple[int, int]:
        """Pone los nombres conocidos (HbS...). Devuelve (encontrados, de los genes guardados)."""
        results: dict[tuple[str, str], bool] = {}      # (gen, nombre) -> encontrado
        with self.conn:
            for symbol, hgvs_c, hgvs_p, name in COMMON_NAMES:
                gene_id = self._gene_id(symbol)
                if gene_id is None:
                    continue
                changed = self.conn.execute("UPDATE variants SET common_name = ? WHERE gene_id = ? "
                                            "AND hgvs_c = ?", (name, gene_id, hgvs_c)).rowcount
                if not changed:
                    changed = self.conn.execute("UPDATE variants SET common_name = ? WHERE gene_id = ? "
                                                "AND hgvs_p = ?", (name, gene_id, hgvs_p)).rowcount
                # Un nombre puede tener varias formas (G6PD): basta con encontrar una
                results[(symbol, name)] = results.get((symbol, name), False) or bool(changed)
        if report is not None:
            report.missing_names.extend(key for key, ok in results.items() if not ok)
        return sum(results.values()), len(results)

    # ------------------------------------------------------------------ todo
    def run(self, genes=CATALOG_GENES, redo: bool = False) -> CatalogReport:
        report = CatalogReport()
        report.removed_genes = self.clean()
        report.repaired = list(self.repaired)
        with self.conn:
            self.conn.execute("INSERT OR REPLACE INTO catalog_info VALUES ('dnainvent_version', ?)",
                              (__version__,))
        for symbol in genes:
            symbol = symbol.upper()
            self.on_gene(symbol)
            actions = {"gene": lambda: self._import_gene(symbol),
                       "variants": lambda: self._variants(symbol),
                       "complete": lambda: self._complete(symbol),
                       "sequences": lambda: self._sequences(symbol),
                       "gene_drugs": lambda: self._gene_drugs(symbol),
                       "disease_drugs": lambda: self._disease_drugs(symbol)}
            for step in GENE_STEPS:
                ok = self._step(report, symbol, step, actions[step], redo)
                if not ok and step in ("gene", "variants"):
                    break        # sin el gen o sus variantes, los demás pasos no tienen sentido
        self.on_gene("*")
        self._step(report, "*", "combined", self._remove_combined, redo)
        # Orphanet y los nombres se rehacen siempre: son rápidos (caché) y dependen de lo nuevo
        self._step(report, "*", "orphanet", self._orphanet, redo=True)
        self._step(report, "*", "common_names",
                   lambda: "{} de {}".format(*self.apply_common_names(report)), redo=True)
        with self.conn:
            self.conn.execute("INSERT OR REPLACE INTO catalog_info VALUES ('built_at', ?)", (now_iso(),))
            self.conn.execute("INSERT OR REPLACE INTO catalog_info VALUES ('genes', ?)",
                              (",".join(catalog_genes(self.conn)),))
        return report
