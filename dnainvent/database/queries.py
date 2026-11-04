"""Consultas de lectura sobre la base de datos: lo que enseñan la consola y (más adelante) la web.

Todas devuelven filas de sqlite3 (se leen como diccionarios: row["symbol"]) o listas de ellas.
"""

import re
import sqlite3

SIGNIFICANCE_ORDER = ("pathogenic", "pathogenic/likely pathogenic", "likely pathogenic",
                      "uncertain significance", "conflicting classifications of pathogenicity",
                      "likely benign", "benign/likely benign", "benign")


def genes(conn: sqlite3.Connection) -> list[sqlite3.Row]:
    """Genes guardados, con cuántas variantes tiene cada uno."""
    return conn.execute(
        """SELECT g.*, o.scientific_name, COUNT(v.id) AS variant_count
           FROM genes g JOIN organisms o ON o.taxon_id = g.taxon_id
           LEFT JOIN variants v ON v.gene_id = g.id
           GROUP BY g.id ORDER BY g.symbol""").fetchall()


def find_gene(conn: sqlite3.Connection, symbol: str, taxon_id: int | None = None) -> sqlite3.Row | None:
    sql = ("SELECT g.*, o.scientific_name FROM genes g JOIN organisms o ON o.taxon_id = g.taxon_id "
           "WHERE g.symbol = ? COLLATE NOCASE")
    params: tuple = (symbol.strip(),)
    if taxon_id is not None:
        sql += " AND g.taxon_id = ?"
        params += (taxon_id,)
    return conn.execute(sql + " ORDER BY g.taxon_id = 9606 DESC", params).fetchone()


def sequences(conn: sqlite3.Connection, gene_id: int) -> list[sqlite3.Row]:
    return conn.execute(
        "SELECT * FROM sequences WHERE gene_id = ? ORDER BY is_reference DESC, mol_type, accession, version",
        (gene_id,)).fetchall()


def reference_mrna(conn: sqlite3.Connection, gene_id: int) -> sqlite3.Row | None:
    return conn.execute("SELECT * FROM sequences WHERE gene_id = ? AND mol_type = 'mrna' "
                        "AND is_reference = 1", (gene_id,)).fetchone()


def significance_counts(conn: sqlite3.Connection, gene_id: int) -> list[tuple[str | None, int]]:
    """(clasificación, nº de variantes), de patogénica a benigna; el resto al final."""
    rows = conn.execute("SELECT clinical_significance AS sig, COUNT(*) AS n FROM variants "
                        "WHERE gene_id = ? GROUP BY sig", (gene_id,)).fetchall()
    rank = {s: i for i, s in enumerate(SIGNIFICANCE_ORDER)}
    return sorted(((r["sig"], r["n"]) for r in rows),
                  key=lambda item: (rank.get(item[0], len(rank)), -item[1]))


def diseases(conn: sqlite3.Connection, gene_id: int, limit: int = 10) -> list[sqlite3.Row]:
    """Enfermedades del gen con el nº de variantes de cada efecto, PARA ESA ENFERMEDAD.
    Primero las que tienen más variantes causales."""
    return conn.execute(
        """SELECT d.id, d.name_en, d.mondo_id,
                  COUNT(CASE WHEN vd.effect = 'causal' THEN 1 END) AS causal,
                  COUNT(CASE WHEN vd.effect = 'risk_factor' THEN 1 END) AS risk_factor,
                  COUNT(CASE WHEN vd.effect = 'protective' THEN 1 END) AS protective,
                  COUNT(*) AS total
           FROM variant_diseases vd JOIN variants v ON v.id = vd.variant_id
           JOIN diseases d ON d.id = vd.disease_id
           WHERE v.gene_id = ?
           GROUP BY d.id
           HAVING causal + risk_factor + protective > 0
           ORDER BY causal DESC, risk_factor + protective DESC, d.name_en
           LIMIT ?""", (gene_id, limit)).fetchall()


def find_variants(conn: sqlite3.Connection, gene_id: int, text: str = "",
                  causal_only: bool = False, limit: int = 20) -> list[sqlite3.Row]:
    """Variantes del gen que coinciden con el texto: c.20A>T, p.Glu7Val, Glu7Val, rs334 o HbS.
    Sin texto, las primeras por posición en el gen."""
    text = text.strip()
    where, params = ["v.gene_id = ?"], [gene_id]
    if text:
        if not text.lower().startswith(("c.", "p.", "rs")) and re.match(r"[A-Za-z]{3}\d", text):
            text = "p." + text
        where.append("(v.hgvs_c = ? OR v.hgvs_p = ? COLLATE NOCASE OR v.dbsnp_rs = ? COLLATE NOCASE "
                     "OR v.common_name = ? COLLATE NOCASE OR v.hgvs_c LIKE ?)")
        params += [text, text, text, text, f"{text}%"]
    if causal_only:
        where.append("EXISTS (SELECT 1 FROM variant_diseases vd WHERE vd.variant_id = v.id "
                     "AND vd.effect = 'causal')")
    return conn.execute(
        f"""SELECT v.*, s.accession, s.version FROM variants v
            JOIN sequences s ON s.id = v.sequence_id
            WHERE {' AND '.join(where)}
            ORDER BY CAST(SUBSTR(v.hgvs_c, 3) AS INTEGER), v.hgvs_c LIMIT ?""",
        params + [limit]).fetchall()


def variant_diseases(conn: sqlite3.Connection, variant_id: int) -> list[sqlite3.Row]:
    """Enfermedades de una variante con la clasificación de cada una."""
    return conn.execute(
        """SELECT d.name_en, d.mondo_id, vd.effect, vd.clinical_significance, vd.review_status,
                  vd.rcv_accession
           FROM variant_diseases vd JOIN diseases d ON d.id = vd.disease_id
           WHERE vd.variant_id = ?
           ORDER BY vd.effect IS NULL, vd.effect, d.name_en""", (variant_id,)).fetchall()


def search(conn: sqlite3.Connection, text: str, limit: int = 20) -> list[sqlite3.Row]:
    """Busca genes y enfermedades por palabras (sin importar acentos ni mayúsculas).
    'anemia' encuentra 'Sickle cell anemia'; 'hemo' encuentra 'hemoglobin' (prefijo)."""
    words = re.findall(r"\w+", text, flags=re.UNICODE)
    if not words:
        return []
    query = " ".join(f'"{w}"*' for w in words)      # cada palabra como prefijo, todas a la vez
    return conn.execute(
        """SELECT si.kind, si.ref_id,
                  COALESCE(g.symbol || ' · ' || g.name, d.name_en) AS label,
                  (SELECT COUNT(*) FROM variant_diseases vd WHERE vd.disease_id = d.id) AS variants
           FROM search_index si
           LEFT JOIN genes g ON si.kind = 'gene' AND g.id = si.ref_id
           LEFT JOIN diseases d ON si.kind = 'disease' AND d.id = si.ref_id
           WHERE search_index MATCH ?
           ORDER BY si.kind DESC, si.rank LIMIT ?""", (query, limit)).fetchall()


def genes_for_disease(conn: sqlite3.Connection, disease_id: int) -> list[sqlite3.Row]:
    """Genes con variantes de esa enfermedad y cuántas causales tiene cada uno."""
    return conn.execute(
        """SELECT g.symbol, g.name, COUNT(CASE WHEN vd.effect = 'causal' THEN 1 END) AS causal,
                  COUNT(*) AS total
           FROM variant_diseases vd JOIN variants v ON v.id = vd.variant_id
           JOIN genes g ON g.id = v.gene_id
           WHERE vd.disease_id = ? GROUP BY g.id ORDER BY causal DESC""", (disease_id,)).fetchall()
