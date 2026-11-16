"""Los genes de serie: la base de datos que DNAinvent trae ya hecha (consola y web).

    builder = CatalogBuilder(open_database("catalog.db"))
    report = builder.run()                         # los 14 genes; se puede cortar y seguir

Se construye en el ordenador de quien publica (tools/build_catalog.py), porque hace miles de
consultas al NCBI y a ChEMBL. La web solo la LEE: así el servidor no llama a ninguna API.

Pasos de cada gen, en orden:
    gene           NCBI Gene: el gen, su ARNm de referencia y su proteína
    variants       ClinVar: todas sus variantes, con la clasificación por enfermedad y la
                   frecuencia en la población
    sequences      las secuencias de los ARNm sobre los que se describen las variantes (para
                   poder comparar sana frente a variante sin conexión)
    gene_drugs     ChEMBL: fármacos que actúan sobre su proteína
    disease_drugs  ChEMBL: fármacos para sus enfermedades principales
Y al final, para todo el catálogo:
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
from .clinvar import ClinVarProvider
from .ncbi import HUMAN, NCBIProvider
from .orphanet import OrphanetProvider
from .records import now_iso

Progress = Callable[[str, dict], None]

CATALOG_GENES = ("HBB", "CFTR", "BRCA1", "BRCA2", "TP53", "HTT", "DMD", "F8", "PAH", "APOE",
                 "LDLR", "SMN1", "G6PD", "LCT")
GENE_STEPS = ("gene", "variants", "sequences", "gene_drugs", "disease_drugs")
DISEASES_WITH_DRUGS = 5       # enfermedades de cada gen para las que se buscan fármacos

# Nombres conocidos de variantes: (gen, HGVS c., HGVS p., nombre). Se buscan primero por c. y,
# si no está, por p. Los nombres son los habituales en la literatura; su explicación en cada
# idioma va aparte (i18n), no aquí.
#   HbS / HbC: hemoglobina de la anemia falciforme y la hemoglobina C (codón 6 en la numeración
#              antigua, sin la metionina inicial; 7 en HGVS).
#   185delAG, 5382insC (BRCA1) y 6174delT (BRCA2): las tres variantes fundadoras asquenazíes,
#              con su nombre antiguo (numeración de BIC).
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

    def _import_gene(self, symbol: str) -> str:
        gene_id = self._gene_id(symbol) or self.ncbi.import_gene(self.conn, symbol, HUMAN, self.progress)
        row = self.conn.execute("SELECT symbol, name FROM genes WHERE id = ?", (gene_id,)).fetchone()
        if (row["symbol"] or "").upper() != symbol.upper():
            found = row["symbol"]
            self.forget_gene(gene_id)
            raise StepIncomplete(f"el NCBI ha devuelto {found} en lugar de {symbol}")
        return row["name"] or ""

    def _variants(self, symbol: str) -> str:
        stats = self.clinvar.import_variants(self.conn, symbol, max_variants=self.max_variants,
                                             progress=self.progress)
        detail = f"{stats.saved} guardadas de {stats.found}; {stats.with_frequency} con frecuencia"
        if stats.condition_batches_failed:
            # Sin la clasificación por enfermedad el gen queda cojo: que se repita
            raise StepIncomplete(f"{stats.condition_batches_failed} lotes de RCV fallidos")
        return detail

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
        found = total = 0
        with self.conn:
            for symbol, hgvs_c, hgvs_p, name in COMMON_NAMES:
                gene_id = self._gene_id(symbol)
                if gene_id is None:
                    continue
                total += 1
                changed = self.conn.execute("UPDATE variants SET common_name = ? WHERE gene_id = ? "
                                            "AND hgvs_c = ?", (name, gene_id, hgvs_c)).rowcount
                if not changed:
                    changed = self.conn.execute("UPDATE variants SET common_name = ? WHERE gene_id = ? "
                                                "AND hgvs_p = ?", (name, gene_id, hgvs_p)).rowcount
                if changed:
                    found += 1
                elif report is not None:
                    report.missing_names.append((symbol, name))
        return found, total

    # ------------------------------------------------------------------ todo
    def run(self, genes=CATALOG_GENES, redo: bool = False) -> CatalogReport:
        report = CatalogReport()
        report.removed_genes = self.clean()
        with self.conn:
            self.conn.execute("INSERT OR REPLACE INTO catalog_info VALUES ('dnainvent_version', ?)",
                              (__version__,))
        for symbol in genes:
            symbol = symbol.upper()
            self.on_gene(symbol)
            actions = {"gene": lambda: self._import_gene(symbol),
                       "variants": lambda: self._variants(symbol),
                       "sequences": lambda: self._sequences(symbol),
                       "gene_drugs": lambda: self._gene_drugs(symbol),
                       "disease_drugs": lambda: self._disease_drugs(symbol)}
            for step in GENE_STEPS:
                ok = self._step(report, symbol, step, actions[step], redo)
                if not ok and step in ("gene", "variants"):
                    break        # sin el gen o sus variantes, los demás pasos no tienen sentido
        self.on_gene("*")
        # Orphanet y los nombres se rehacen siempre: son rápidos (caché) y dependen de lo nuevo
        self._step(report, "*", "orphanet", self._orphanet, redo=True)
        self._step(report, "*", "common_names",
                   lambda: "{} de {}".format(*self.apply_common_names(report)), redo=True)
        with self.conn:
            self.conn.execute("INSERT OR REPLACE INTO catalog_info VALUES ('built_at', ?)", (now_iso(),))
            self.conn.execute("INSERT OR REPLACE INTO catalog_info VALUES ('genes', ?)",
                              (",".join(catalog_genes(self.conn)),))
        return report
