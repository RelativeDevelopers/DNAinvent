"""Base de datos local de genes, variantes y enfermedades reales (SQLite, sin dependencias).

    conn = open_database()                       # ~/.dnainvent/dnainvent.db
    NCBIProvider().import_gene(conn, "HBB")      # gen + ARNm + proteína   (ncbi.py)
    ClinVarProvider().import_variants(conn, "HBB")   # variantes y enfermedades (clinvar.py)

Idea general (ver schema.sql):
  * En local se guardan los METADATOS: ocupan poco y sobre ellos se filtra.
  * La secuencia se descarga BAJO DEMANDA (NCBIProvider.ensure_residues) y se guarda.
  * Cada fila sabe de qué fuente viene y con qué licencia (tabla sources).

Como todo el paquete dnainvent, no imprime nada: avisa con progress(clave, datos), con claves
de i18n (database.stage.*), y los fallos son errores con código.

Para usar otro archivo: variable de entorno DNAINVENT_DB. Borrarlo es seguro: se vuelve a
crear vacío y se reimporta lo que haga falta.
"""

import os
import sqlite3
from pathlib import Path

SCHEMA_PATH = Path(__file__).with_name("schema.sql")
SCHEMA_VERSION = 2

# Fuentes conocidas: código -> (nombre, url, licencia, uso comercial permitido)
KNOWN_SOURCES = {
    "refseq": ("NCBI RefSeq", "https://www.ncbi.nlm.nih.gov/refseq/", "public-domain", 1),
    "ncbi_gene": ("NCBI Gene", "https://www.ncbi.nlm.nih.gov/gene/", "public-domain", 1),
    "clinvar": ("NCBI ClinVar", "https://www.ncbi.nlm.nih.gov/clinvar/", "public-domain", 1),
    "mondo": ("Mondo Disease Ontology", "https://mondo.monarchinitiative.org/", "CC-BY-4.0", 1),
    "chembl": ("ChEMBL", "https://www.ebi.ac.uk/chembl/", "CC-BY-SA-3.0", 1),
    "manual": ("DNAinvent (añadido a mano)", None, "public-domain", 1),
}

# Columnas añadidas en cada versión del esquema (para actualizar bases de datos antiguas)
_MIGRATIONS = {
    2: [
        "ALTER TABLE variants ADD COLUMN review_status TEXT",
        "ALTER TABLE variant_diseases ADD COLUMN clinical_significance TEXT",
        "ALTER TABLE variant_diseases ADD COLUMN review_status TEXT",
        "ALTER TABLE variant_diseases ADD COLUMN rcv_accession TEXT",
    ],
}


def default_path() -> Path:
    return Path(os.environ.get("DNAINVENT_DB") or Path.home() / ".dnainvent" / "dnainvent.db")


def open_database(path: str | Path | None = None) -> sqlite3.Connection:
    """Abre (o crea) la base de datos y la deja en la versión actual del esquema."""
    path = default_path() if path is None else path
    if str(path) != ":memory:":
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(path)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    exists = conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'genes'").fetchone()
    if not exists:
        conn.executescript(SCHEMA_PATH.read_text(encoding="utf-8"))
    else:
        _migrate(conn)
    return conn


def _migrate(conn: sqlite3.Connection) -> None:
    version = conn.execute("PRAGMA user_version").fetchone()[0] or 1
    for target in range(version + 1, SCHEMA_VERSION + 1):
        with conn:
            for statement in _MIGRATIONS.get(target, []):
                try:
                    conn.execute(statement)
                except sqlite3.OperationalError as error:
                    if "duplicate column" not in str(error):
                        raise
            conn.execute(f"PRAGMA user_version = {target}")


def get_or_create_source(conn: sqlite3.Connection, code: str) -> int:
    """Id de una fuente de la tabla sources, creándola si aún no existe."""
    row = conn.execute("SELECT id FROM sources WHERE code = ?", (code,)).fetchone()
    if row:
        return row["id"]
    if code not in KNOWN_SOURCES:
        raise ValueError(f"Fuente desconocida: {code}")
    name, url, license_, commercial_ok = KNOWN_SOURCES[code]
    return conn.execute(
        "INSERT INTO sources (code, name, url, license, commercial_ok) VALUES (?, ?, ?, ?, ?)",
        (code, name, url, license_, commercial_ok),
    ).lastrowid
