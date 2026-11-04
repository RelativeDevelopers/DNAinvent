"""NCBI Gene + RefSeq: un gen, su ARNm de referencia y su proteína.

    provider = NCBIProvider()
    gene_id = provider.import_gene(conn, "HBB")            # humano (taxón 9606) por defecto
    gene_id = provider.import_gene(conn, "Hbb-b1", 10090)  # ratón

1. Gene      esearch + esummary: nombre, cromosoma y coordenadas del gen.
2. ARNm      el RefSeq Select (la versión de referencia del gen); si no hay, el primer NM_.
3. Registro  efetch del ARNm en XML INSDSeq (rettype=gbc): la secuencia, dónde empieza y
             acaba el CDS y la proteína que produce (NP_), ya traducida por el NCBI.

Coordenadas: el NCBI las da desde 0 y, si el gen está en la cadena -, al revés (start > stop);
aquí se guardan desde 1, con start <= end y la cadena aparte.
"""

import re
import sqlite3
import xml.etree.ElementTree as ET
from typing import Callable

from ..errors import ServiceError
from . import get_or_create_source
from .eutils import EUtilsClient
from .records import now_iso

HUMAN = 9606
Progress = Callable[[str, dict], None]

# Tipo de molécula de GenBank -> sequences.mol_type
MOLTYPE_MAP = {"mrna": "mrna", "dna": "genomic", "rna": "ncrna", "ncrna": "ncrna",
               "rrna": "ncrna", "trna": "ncrna", "aa": "protein"}


def _no_progress(key: str, data: dict) -> None:
    pass


def parse_location(location: str) -> tuple[int, int] | None:
    """'51..494' -> (51, 494). None si es compuesta (join, complement...): en un ARNm de
    RefSeq el CDS siempre es un tramo continuo."""
    match = re.fullmatch(r"<?(\d+)\.\.>?(\d+)", location.strip())
    return (int(match.group(1)), int(match.group(2))) if match else None


def parse_insdseq(raw: bytes, accession: str, service: str = "NCBI") -> dict:
    """Un registro de nuccore en XML. Entiende INSDSeq (rettype=gbc) y GBSeq (rettype=gb),
    que tienen el mismo contenido con otro prefijo en las etiquetas."""
    try:
        root = ET.fromstring(raw)
    except ET.ParseError:
        raise ServiceError("service_unavailable", service=service) from None
    for p in ("INSD", "GB"):
        seq = root if root.tag == f"{p}Seq" else root.find(f"{p}Seq")
        if seq is not None:
            break
    else:
        raise ServiceError("sequence_not_found", accession=accession)

    acc, _, version = seq.findtext(f"{p}Seq_accession-version", "").partition(".")
    record = {
        "accession": acc,
        "version": int(version) if version else None,
        "mol_type": MOLTYPE_MAP.get(seq.findtext(f"{p}Seq_moltype", "").lower(), "genomic"),
        "description": seq.findtext(f"{p}Seq_definition"),
        "residues": seq.findtext(f"{p}Seq_sequence", "").upper(),
        "taxon_id": None, "cds_start": None, "cds_end": None,
        "genetic_code": None, "protein": None,
    }
    record["length"] = len(record["residues"])

    for feature in seq.iter(f"{p}Feature"):
        key = feature.findtext(f"{p}Feature_key")
        quals: dict[str, list[str]] = {}
        for q in feature.iter(f"{p}Qualifier"):
            quals.setdefault(q.findtext(f"{p}Qualifier_name"), []).append(
                q.findtext(f"{p}Qualifier_value") or "")

        if key == "source":
            for xref in quals.get("db_xref", []):
                if xref.startswith("taxon:"):
                    record["taxon_id"] = int(xref.split(":", 1)[1])
        elif key == "CDS" and record["cds_start"] is None:
            span = parse_location(feature.findtext(f"{p}Feature_location", ""))
            if span:
                record["cds_start"], record["cds_end"] = span
            if "transl_table" in quals:
                record["genetic_code"] = int(quals["transl_table"][0])
            if "protein_id" in quals and "translation" in quals:
                p_acc, _, p_ver = quals["protein_id"][0].partition(".")
                translation = quals["translation"][0].replace(" ", "").upper()
                record["protein"] = {
                    "accession": p_acc, "version": int(p_ver) if p_ver else None,
                    "residues": translation, "length": len(translation),
                    "description": (quals.get("product") or [None])[0],
                }
    return record


class NCBIProvider(EUtilsClient):
    service = "NCBI"

    # ------------------------------------------------------------------ Gene
    def search_gene(self, symbol: str, taxon_id: int = HUMAN) -> int | None:
        data = self.json("esearch.fcgi", {"db": "gene", "term": f"{symbol}[sym] AND {taxon_id}[taxid]"})
        ids = data.get("esearchresult", {}).get("idlist", [])
        return int(ids[0]) if ids else None

    def fetch_gene_summary(self, gene_id: int) -> dict:
        data = self.json("esummary.fcgi", {"db": "gene", "id": gene_id})
        doc = data.get("result", {}).get(str(gene_id))
        if not doc:
            raise ServiceError("search_failed", service=self.service)
        organism = doc.get("organism", {})
        summary = {
            "ncbi_gene_id": int(gene_id), "symbol": doc.get("name"),
            "name": doc.get("description"), "chromosome": doc.get("chromosome") or None,
            "taxon_id": int(organism.get("taxid", 0)) or None,
            "scientific_name": organism.get("scientificname"),
            "common_name": organism.get("commonname"),
            "start_pos": None, "end_pos": None, "strand": None, "assembly": None,
        }
        genomic = doc.get("genomicinfo") or []
        if genomic:
            info = genomic[0]
            start, stop = int(info["chrstart"]), int(info["chrstop"])
            summary["start_pos"] = min(start, stop) + 1
            summary["end_pos"] = max(start, stop) + 1
            summary["strand"] = -1 if start > stop else 1
            summary["assembly"] = info.get("chraccver")     # p. ej. NC_000011.10
        return summary

    def find_reference_transcript(self, symbol: str, taxon_id: int, gene_id: int) -> str | None:
        """Accesión.versión del ARNm de referencia: RefSeq Select o, si no hay, el primer NM_."""
        term = (f"{symbol}[Gene Name] AND txid{taxon_id}[Organism:noexp] "
                "AND refseq_select[filter] AND biomol_mrna[prop]")
        uids = self.json("esearch.fcgi", {"db": "nuccore", "term": term}) \
            .get("esearchresult", {}).get("idlist", [])
        if not uids:
            links = self.json("elink.fcgi", {"dbfrom": "gene", "db": "nuccore", "id": gene_id,
                                             "linkname": "gene_nuccore_refseqrna"})
            for linkset in links.get("linksets", []):
                for db in linkset.get("linksetdbs", []):
                    uids.extend(db.get("links", []))
        if not uids:
            return None
        result = self.json("esummary.fcgi", {"db": "nuccore", "id": ",".join(map(str, uids))}) \
            .get("result", {})
        accessions = [result[str(uid)]["accessionversion"] for uid in result.get("uids", [])
                      if "accessionversion" in result.get(str(uid), {})]
        mrnas = sorted(a for a in accessions if a.startswith("NM_"))
        return mrnas[0] if mrnas else (accessions[0] if accessions else None)

    # ------------------------------------------------------------------ Secuencias
    def fetch_sequence(self, accession_version: str) -> dict:
        raw = self.request("efetch.fcgi", {"db": "nuccore", "id": accession_version,
                                           "rettype": "gbc", "retmode": "xml"})
        return parse_insdseq(raw, accession_version, self.service)

    # ------------------------------------------------------------------ Guardar
    def import_gene(self, conn: sqlite3.Connection, symbol: str, taxon_id: int = HUMAN,
                    progress: Progress = _no_progress) -> int:
        """Importa (o actualiza) un gen con su ARNm de referencia y su proteína.
        Devuelve el id local del gen."""
        progress("database.stage.gene", {"symbol": symbol})
        gene_id = self.search_gene(symbol, taxon_id)
        if gene_id is None:
            raise ServiceError("gene_not_found", symbol=symbol, taxon=taxon_id)
        summary = self.fetch_gene_summary(gene_id)
        transcript = self.find_reference_transcript(summary["symbol"], taxon_id, gene_id)
        mrna = None
        if transcript:
            progress("database.stage.sequence", {"accession": transcript})
            mrna = self.fetch_sequence(transcript)

        with conn:          # una sola transacción: o se guarda todo o nada
            gene_source = get_or_create_source(conn, "ncbi_gene")
            seq_source = get_or_create_source(conn, "refseq")
            _upsert_organism(conn, summary)
            local_gene_id = _upsert_gene(conn, summary, gene_source)
            if mrna:
                mrna_id = upsert_sequence(conn, mrna, local_gene_id, summary["taxon_id"],
                                          seq_source, is_reference=1)
                if mrna["protein"]:
                    protein = dict(mrna["protein"], mol_type="protein", cds_start=None,
                                   cds_end=None, genetic_code=mrna["genetic_code"])
                    protein_id = upsert_sequence(conn, protein, local_gene_id,
                                                 summary["taxon_id"], seq_source, is_reference=1)
                    conn.execute("INSERT OR IGNORE INTO sequence_links VALUES (?, ?, 'translates_to')",
                                 (mrna_id, protein_id))
            conn.execute("DELETE FROM search_index WHERE kind = 'gene' AND ref_id = ?",
                         (local_gene_id,))
            conn.execute("INSERT INTO search_index (kind, ref_id, text) VALUES ('gene', ?, ?)",
                         (local_gene_id, f"{summary['symbol']} {summary['name'] or ''}"))
        return local_gene_id

    def ensure_residues(self, conn: sqlite3.Connection, sequence_id: int) -> str:
        """Secuencia bajo demanda: si solo están los metadatos, la descarga y la guarda.
        Sirve también para versiones antiguas de un ARNm (NM_000518.4) que citan las variantes."""
        row = conn.execute("SELECT * FROM sequences WHERE id = ?", (sequence_id,)).fetchone()
        if row is None:
            raise ServiceError("sequence_not_found", accession=str(sequence_id))
        if row["residues"]:
            return row["residues"]
        accession = f"{row['accession']}.{row['version']}"
        if row["mol_type"] == "protein":        # las proteínas llegan junto a su ARNm
            raise ServiceError("sequence_not_found", accession=accession)
        record = self.fetch_sequence(accession)
        with conn:
            conn.execute(
                """UPDATE sequences SET residues = ?, length = ?, cds_start = ?, cds_end = ?,
                   genetic_code = COALESCE(?, genetic_code), description = COALESCE(?, description),
                   fetched_at = ? WHERE id = ?""",
                (record["residues"], record["length"], record["cds_start"], record["cds_end"],
                 record["genetic_code"], record["description"], now_iso(), sequence_id))
        return record["residues"]


# ---------------------------------------------------------------------- escritura
def _upsert_organism(conn, summary):
    conn.execute(
        """INSERT INTO organisms (taxon_id, rank, scientific_name, name_en)
           VALUES (?, 'species', ?, ?) ON CONFLICT (taxon_id) DO NOTHING""",
        (summary["taxon_id"], summary["scientific_name"] or "?", summary["common_name"]))


def _upsert_gene(conn, summary, source_id) -> int:
    row = conn.execute("SELECT id FROM genes WHERE ncbi_gene_id = ?",
                       (summary["ncbi_gene_id"],)).fetchone()
    values = (summary["taxon_id"], summary["symbol"], summary["name"], summary["chromosome"],
              summary["start_pos"], summary["end_pos"], summary["strand"], summary["assembly"],
              source_id)
    if row:
        conn.execute(
            """UPDATE genes SET taxon_id = ?, symbol = ?, name = ?, chromosome = ?, start_pos = ?,
               end_pos = ?, strand = ?, assembly = ?, source_id = ? WHERE id = ?""",
            values + (row["id"],))
        return row["id"]
    return conn.execute(
        """INSERT INTO genes (taxon_id, symbol, name, chromosome, start_pos, end_pos, strand,
                              assembly, source_id, ncbi_gene_id)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        values + (summary["ncbi_gene_id"],)).lastrowid


def upsert_sequence(conn, record, gene_id, taxon_id, source_id, is_reference=0) -> int:
    row = conn.execute("SELECT id FROM sequences WHERE accession = ? AND version = ?",
                       (record["accession"], record["version"])).fetchone()
    values = (gene_id, taxon_id, record["mol_type"], is_reference, record["cds_start"],
              record["cds_end"], record["genetic_code"], record["length"],
              record["description"], record["residues"], now_iso(), source_id)
    if row:
        conn.execute(
            """UPDATE sequences SET gene_id = ?, taxon_id = ?, mol_type = ?, is_reference = ?,
               cds_start = ?, cds_end = ?, genetic_code = ?, length = ?, description = ?,
               residues = ?, fetched_at = ?, source_id = ? WHERE id = ?""",
            values + (row["id"],))
        return row["id"]
    return conn.execute(
        """INSERT INTO sequences (gene_id, taxon_id, mol_type, is_reference, cds_start, cds_end,
                                  genetic_code, length, description, residues, fetched_at,
                                  source_id, accession, version)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        values + (record["accession"], record["version"])).lastrowid
