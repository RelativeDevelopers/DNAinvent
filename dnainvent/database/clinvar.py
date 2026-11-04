"""ClinVar: las variantes de un gen humano, su clasificación y sus enfermedades.

    stats = ClinVarProvider().import_variants(conn, "HBB")

1. Buscar     esearch "HBB[gene]": todas las variantes asociadas al gen. No se usa
              single_gene[prop]: ClinVar asocia la mayoría también a regiones reguladoras
              vecinas (LOC...) y ese filtro las descarta casi todas (26 de 1.974 en HBB).
2. Resúmenes  esummary JSON, de 200 en 200: el nombre HGVS (NM_000518.5(HBB):c.20A>T
              (p.Glu7Val)), la clasificación GLOBAL de la variante y las enfermedades
              estudiadas, con sus IDs de MedGen y MONDO.
3. Por enfermedad  efetch rettype=vcv (XML), de 25 en 25: cada registro RCV es una pareja
              variante + enfermedad con SU PROPIA clasificación.

Por qué hace falta el paso 3: la clasificación global junta todo lo que han enviado los
laboratorios. HbS es "patogénica" en global, pero para la malaria es "protectora"; si se usara
la global para cada enfermedad, saldría que HbS causa la malaria. Así que el efecto
(causal / factor de riesgo / protectora) se saca de cada RCV. Si el RCV no se puede leer, la
enfermedad se enlaza igualmente, pero sin efecto (NULL): se sabe que se estudió, no qué hace.

Se guardan solo las variantes descritas sobre un ARNm de este gen (las que DNAinvent puede
situar en la secuencia). Las grandes deleciones de cromosoma y las de genes vecinos se
cuentan, pero no se guardan.
"""

import re
import sqlite3
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import Callable

from ..errors import ServiceError
from . import get_or_create_source
from .eutils import EUtilsClient
from .ncbi import HUMAN, NCBIProvider

Progress = Callable[[str, dict], None]
SEARCH_PAGE = 5000        # IDs por petición esearch
SUMMARY_BATCH = 200       # IDs por petición esummary
VCV_BATCH = 25            # IDs por petición efetch (el XML de algunas variantes es grande)

# NM_000518.5(HBB):c.20A>T (p.Glu7Val)
HGVS_NAME = re.compile(r"(?P<acc>[NX][MR]_\d+)\.(?P<ver>\d+)\((?P<gene>[^)]+)\):"
                       r"(?P<c>[cn]\.\S+?)(?:\s+\((?P<p>p\.[^)]+)\))?\s*$")
SIMPLE_SNV = re.compile(r"c\.-?\*?\d+([ACGT])>([ACGT])")

VARIANT_TYPES = {"single nucleotide variant": "snv", "deletion": "deletion",
                 "insertion": "insertion", "indel": "indel", "duplication": "duplication",
                 "copy number loss": "cnv", "copy number gain": "cnv"}

# Rasgos que ClinVar usa cuando no hay una enfermedad concreta
NON_DISEASE_TRAITS = {"not provided", "not specified", "see cases",
                      "inborn genetic diseases", ""}

EFFECT_STRENGTH = {"causal": 3, "risk_factor": 2, "protective": 1, None: 0}


def _no_progress(key: str, data: dict) -> None:
    pass


@dataclass
class ImportStats:
    found: int = 0               # variantes que ClinVar asocia al gen
    saved: int = 0
    pathogenic: int = 0          # patogénicas o probablemente patogénicas (global)
    not_on_mrna: int = 0         # grandes deleciones, cromosoma... (no se guardan)
    other_gene: int = 0          # descritas sobre el ARNm de un gen vecino
    filtered: int = 0            # descartadas por only_pathogenic
    with_conditions: int = 0     # variantes con clasificación por enfermedad (RCV)
    condition_batches_failed: int = 0   # lotes de RCV que no llegaron o no se pudieron leer
    unexpected_records: int = 0  # registros que ClinVar devolvió sin haberlos pedido (descartados)
    new_diseases: int = 0
    gene_imported: bool = False


@dataclass(frozen=True)
class Condition:
    """Un registro RCV: esta variante, para esta enfermedad, tiene esta clasificación."""

    rcv: str
    name: str
    medgen_id: str | None
    significance: str | None
    review_status: str | None


# ---------------------------------------------------------------------- lectura
def parse_variation_name(name: str | None) -> dict | None:
    """'NM_000518.5(HBB):c.20A>T (p.Glu7Val)' -> dict; None si no está descrita en un ARNm."""
    match = HGVS_NAME.match(name or "")
    if not match:
        return None
    return {"accession": match.group("acc"), "version": int(match.group("ver")),
            "gene": match.group("gene"), "hgvs_c": match.group("c"), "hgvs_p": match.group("p")}


def classify_effect(significance: str | None) -> str | None:
    """Qué hace la variante en una enfermedad, a partir de la clasificación de ClinVar."""
    s = (significance or "").lower()
    if "risk factor" in s or "risk allele" in s:
        return "risk_factor"
    if "protective" in s:
        return "protective"
    if "conflicting" in s or "benign" in s:
        return None
    if "pathogenic" in s or "oncogenic" in s:      # oncogénica: impulsa un cáncer (somática)
        return "causal"
    return None


def parse_summary(uid: str, doc: dict) -> dict | None:
    """Un documento esummary de ClinVar, listo para guardar (None si no está en un ARNm)."""
    vset = (doc.get("variation_set") or [{}])[0]
    hgvs = parse_variation_name(vset.get("variation_name")) or parse_variation_name(doc.get("title"))
    if hgvs is None:
        return None

    germline = doc.get("germline_classification") or {}
    somatic = ((doc.get("oncogenicity_classification") or {}).get("description")
               or (doc.get("clinical_impact_classification") or {}).get("description"))
    significance = germline.get("description") or somatic or None
    origin = "germline" if germline.get("description") else ("somatic" if somatic else "unknown")

    dbsnp = next((f"rs{x['db_id']}" for x in vset.get("variation_xrefs") or []
                  if x.get("db_source") == "dbSNP" and x.get("db_id")), None)
    traits = []
    for trait in germline.get("trait_set") or []:
        name = (trait.get("trait_name") or "").strip()
        if name.lower() in NON_DISEASE_TRAITS:
            continue
        xrefs = {x.get("db_source"): x.get("db_id") for x in trait.get("trait_xrefs") or []}
        traits.append({"name": name, "mondo_id": xrefs.get("MONDO"),
                       "medgen_id": xrefs.get("MedGen")})

    snv = SIMPLE_SNV.fullmatch(hgvs["hgvs_c"])
    return dict(
        hgvs, clinvar_id=int(uid),
        variant_type=VARIANT_TYPES.get((doc.get("obj_type") or "").lower(), "other"),
        clinical_significance=significance.lower() if significance else None,
        review_status=germline.get("review_status") or None,
        origin=origin, dbsnp_rs=dbsnp,
        ref_allele=snv.group(1) if snv else None, alt_allele=snv.group(2) if snv else None,
        traits=traits,
    )


def parse_vcv(raw: bytes) -> dict[int, list[Condition]]:
    """XML VCV de ClinVar -> {ID de variante: [registros RCV]}.

    <VariationArchive VariationID="15333" ...>
      ... <RCVList>
            <RCVAccession Accession="RCV000016574" ...>
              <ClassifiedConditionList>
                <ClassifiedCondition DB="MedGen" ID="C0002895">Sickle cell disease</...>
              <RCVClassifications><GermlineClassification>
                <ReviewStatus>...</ReviewStatus><Description>Pathogenic</Description>
    """
    try:
        root = ET.fromstring(raw)
    except ET.ParseError:
        raise ServiceError("service_unavailable", service="ClinVar") from None
    archives = [root] if root.tag == "VariationArchive" else root.iter("VariationArchive")
    out: dict[int, list[Condition]] = {}
    for archive in archives:
        try:
            variation_id = int(archive.get("VariationID", ""))
        except ValueError:
            continue
        conditions = out.setdefault(variation_id, [])
        for rcv in archive.iter("RCVAccession"):
            classification = None
            for tag in ("GermlineClassification", "OncogenicityClassification",
                        "SomaticClinicalImpact"):
                classification = rcv.find(f"RCVClassifications/{tag}")
                if classification is not None:
                    break
            significance = review = None
            if classification is not None:
                significance = (classification.findtext("Description") or "").strip() or None
                review = (classification.findtext("ReviewStatus") or "").strip() or None
            for condition in rcv.iter("ClassifiedCondition"):
                name = (condition.text or "").strip()
                medgen = condition.get("ID") if condition.get("DB") == "MedGen" else None
                conditions.append(Condition(rcv.get("Accession", ""), name, medgen,
                                            significance.lower() if significance else None,
                                            review))
    return out


# ---------------------------------------------------------------------- proveedor
class ClinVarProvider(EUtilsClient):
    service = "ClinVar"

    def search_variant_ids(self, symbol: str, max_variants: int | None = None) -> list[str]:
        ids: list[str] = []
        while max_variants is None or len(ids) < max_variants:
            page = SEARCH_PAGE if max_variants is None else min(SEARCH_PAGE, max_variants - len(ids))
            result = self.json("esearch.fcgi", {"db": "clinvar", "term": f"{symbol}[gene]",
                                                "retstart": len(ids), "retmax": page})["esearchresult"]
            batch = result.get("idlist", [])
            ids.extend(batch)
            if not batch or len(ids) >= int(result.get("count", 0)):
                break
        return ids

    def fetch_summaries(self, ids: list[str], progress: Progress = _no_progress) -> list[tuple[str, dict]]:
        docs = []
        for i in range(0, len(ids), SUMMARY_BATCH):
            chunk = ids[i:i + SUMMARY_BATCH]
            result = self.json("esummary.fcgi", {"db": "clinvar", "id": ",".join(chunk)}).get("result", {})
            docs.extend((uid, result[uid]) for uid in result.get("uids", []) if uid in result)
            progress("database.stage.summaries", {"done": i + len(chunk), "total": len(ids)})
        return docs

    def fetch_conditions(self, ids: list[int], progress: Progress = _no_progress,
                         stats: ImportStats | None = None) -> dict[int, list[Condition]]:
        """Registros RCV de cada variante. Si un lote falla, sus variantes quedan sin
        clasificación por enfermedad, pero el resto de la importación sigue (y se cuenta en
        stats, para avisar).

        is_variationid=true es imprescindible: sin él, ClinVar no lee los números como IDs de
        variante y devuelve OTRAS variantes (al pedir unas de TP53 llegaban de MSH2 y BRCA1).
        Por si acaso, solo se aceptan los registros que se han pedido."""
        stats = stats if stats is not None else ImportStats()
        out: dict[int, list[Condition]] = {}
        for i in range(0, len(ids), VCV_BATCH):
            chunk = ids[i:i + VCV_BATCH]
            try:
                raw = self.request("efetch.fcgi", {"db": "clinvar", "rettype": "vcv",
                                                   "is_variationid": "true",
                                                   "id": ",".join(map(str, chunk))})
                records = parse_vcv(raw)
            except ServiceError:
                stats.condition_batches_failed += 1
            else:
                wanted = set(chunk)
                stats.unexpected_records += sum(1 for v in records if v not in wanted)
                out.update((v, c) for v, c in records.items() if v in wanted)
            progress("database.stage.conditions", {"done": i + len(chunk), "total": len(ids)})
        return out

    def import_variants(self, conn: sqlite3.Connection, symbol: str, only_pathogenic: bool = False,
                        max_variants: int | None = None, progress: Progress = _no_progress,
                        gene_provider: NCBIProvider | None = None) -> ImportStats:
        """Importa las variantes de ClinVar de un gen humano (y el gen, si aún no está)."""
        stats = ImportStats()
        gene = conn.execute("SELECT id, symbol FROM genes WHERE symbol = ? COLLATE NOCASE "
                            "AND taxon_id = ?", (symbol, HUMAN)).fetchone()
        if gene is None:
            gene_provider = gene_provider or NCBIProvider(http=self.http, email=self.email,
                                                          api_key=self.api_key, sleep=self.sleep)
            gene_id = gene_provider.import_gene(conn, symbol, HUMAN, progress)
            stats.gene_imported = True
        else:
            gene_id = gene["id"]
        symbol = conn.execute("SELECT symbol FROM genes WHERE id = ?", (gene_id,)).fetchone()[0]

        progress("database.stage.variants", {"symbol": symbol})
        ids = self.search_variant_ids(symbol, max_variants)
        stats.found = len(ids)

        variants = []
        for uid, doc in self.fetch_summaries(ids, progress):
            variant = parse_summary(uid, doc)
            if variant is None:
                stats.not_on_mrna += 1
            elif variant["gene"].upper() != symbol.upper():
                stats.other_gene += 1
            elif only_pathogenic and classify_effect(variant["clinical_significance"]) != "causal":
                stats.filtered += 1
            else:
                variants.append(variant)

        conditions = self.fetch_conditions([v["clinvar_id"] for v in variants], progress, stats)
        progress("database.stage.saving", {"count": len(variants)})
        with conn:
            clinvar_src = get_or_create_source(conn, "clinvar")
            refseq_src = get_or_create_source(conn, "refseq")
            for variant in variants:
                seq_id = _ensure_transcript(conn, variant, gene_id, refseq_src)
                variant_id = _upsert_variant(conn, variant, gene_id, seq_id, clinvar_src)
                stats.saved += 1
                stats.pathogenic += classify_effect(variant["clinical_significance"]) == "causal"
                rcvs = conditions.get(variant["clinvar_id"])
                if rcvs:
                    stats.with_conditions += 1
                stats.new_diseases += _link_diseases(conn, variant_id, gene_id, variant["traits"],
                                                     rcvs, clinvar_src)
            # Quita los enlaces gen-enfermedad que ya no respalda ninguna variante (p. ej. los
            # de una importación anterior que usaba la clasificación global)
            rebuild_gene_diseases(conn, gene_id)
        return stats


# ---------------------------------------------------------------------- escritura
def _ensure_transcript(conn, variant, gene_id, source_id) -> int:
    """El ARNm sobre el que se describe la variante. Si no está, solo sus metadatos: la
    secuencia se descarga cuando haga falta (NCBIProvider.ensure_residues)."""
    row = conn.execute("SELECT id FROM sequences WHERE accession = ? AND version = ?",
                       (variant["accession"], variant["version"])).fetchone()
    if row:
        return row["id"]
    mol_type = "mrna" if variant["accession"][1] == "M" else "ncrna"
    return conn.execute(
        "INSERT INTO sequences (gene_id, taxon_id, accession, version, mol_type, source_id) "
        "VALUES (?, ?, ?, ?, ?, ?)",
        (gene_id, HUMAN, variant["accession"], variant["version"], mol_type, source_id)).lastrowid


def _upsert_variant(conn, v, gene_id, seq_id, source_id) -> int:
    # 1) ya importada de ClinVar; 2) añadida a mano con la misma HGVS: se completa y se
    #    conserva su nombre común (HbS)
    row = conn.execute("SELECT id FROM variants WHERE clinvar_id = ?", (v["clinvar_id"],)).fetchone()
    if row is None:
        row = conn.execute("SELECT id FROM variants WHERE sequence_id = ? AND hgvs_c = ? "
                           "AND clinvar_id IS NULL", (seq_id, v["hgvs_c"])).fetchone()
    values = (gene_id, seq_id, v["hgvs_c"], v["hgvs_p"], v["variant_type"], v["ref_allele"],
              v["alt_allele"], v["origin"], v["clinical_significance"], v["review_status"],
              v["dbsnp_rs"], v["clinvar_id"], source_id)
    if row:
        conn.execute(
            """UPDATE variants SET gene_id = ?, sequence_id = ?, hgvs_c = ?, hgvs_p = ?,
               variant_type = ?, ref_allele = COALESCE(?, ref_allele),
               alt_allele = COALESCE(?, alt_allele), origin = ?, clinical_significance = ?,
               review_status = ?, dbsnp_rs = COALESCE(?, dbsnp_rs), clinvar_id = ?, source_id = ?
               WHERE id = ?""", values + (row["id"],))
        return row["id"]
    return conn.execute(
        """INSERT INTO variants (gene_id, sequence_id, hgvs_c, hgvs_p, variant_type, ref_allele,
                                 alt_allele, origin, clinical_significance, review_status,
                                 dbsnp_rs, clinvar_id, source_id)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""", values).lastrowid


def _link_diseases(conn, variant_id, gene_id, traits, rcvs, source_id) -> int:
    """Enlaza la variante con sus enfermedades. Devuelve cuántas enfermedades son nuevas.

    Con RCV: una fila por enfermedad, con la clasificación de esa enfermedad (si hay varios
    RCV para la misma, gana el efecto más fuerte). Sin RCV: las enfermedades del resumen,
    sin efecto. Los enlaces antiguos de esta variante se rehacen: así una reimportación
    corrige los que se guardaron con la clasificación global."""
    conn.execute("DELETE FROM variant_diseases WHERE variant_id = ?", (variant_id,))
    by_medgen = {t["medgen_id"]: t for t in traits if t["medgen_id"]}
    by_name = {t["name"].lower(): t for t in traits}

    links: dict[int, tuple] = {}
    created = 0
    if rcvs:
        candidates = []
        for c in rcvs:
            if c.name.lower() in NON_DISEASE_TRAITS:
                continue
            trait = by_medgen.get(c.medgen_id) or by_name.get(c.name.lower()) or {}
            candidates.append(({"name": trait.get("name") or c.name,
                                "mondo_id": trait.get("mondo_id")},
                               classify_effect(c.significance), c.significance,
                               c.review_status, c.rcv))
    else:
        candidates = [(t, None, None, None, None) for t in traits]

    for trait, effect, significance, review, rcv in candidates:
        disease_id, new = _upsert_disease(conn, trait)
        created += new
        current = links.get(disease_id)
        if current is None or EFFECT_STRENGTH[effect] > EFFECT_STRENGTH[current[0]]:
            links[disease_id] = (effect, significance, review, rcv)

    for disease_id, (effect, significance, review, rcv) in links.items():
        conn.execute(
            """INSERT INTO variant_diseases (variant_id, disease_id, effect, clinical_significance,
                                             review_status, rcv_accession, source_id)
               VALUES (?, ?, ?, ?, ?, ?, ?)""",
            (variant_id, disease_id, effect, significance, review, rcv, source_id))
        if effect == "causal":
            conn.execute("INSERT OR IGNORE INTO gene_diseases VALUES (?, ?, ?)",
                         (gene_id, disease_id, source_id))
    return created


def _upsert_disease(conn, trait) -> tuple[int, int]:
    """(id, 1 si es nueva) de una enfermedad, buscándola por MONDO o por nombre."""
    row = None
    if trait.get("mondo_id"):
        row = conn.execute("SELECT id FROM diseases WHERE mondo_id = ?", (trait["mondo_id"],)).fetchone()
    if row is None:
        row = conn.execute("SELECT id, mondo_id FROM diseases WHERE name_en = ? COLLATE NOCASE",
                           (trait["name"],)).fetchone()
        if row is not None and trait.get("mondo_id") and not row["mondo_id"]:
            conn.execute("UPDATE diseases SET mondo_id = ? WHERE id = ?", (trait["mondo_id"], row["id"]))
    if row is not None:
        return row["id"], 0
    disease_id = conn.execute("INSERT INTO diseases (mondo_id, name_en) VALUES (?, ?)",
                              (trait.get("mondo_id"), trait["name"])).lastrowid
    conn.execute("INSERT INTO search_index (kind, ref_id, text) VALUES ('disease', ?, ?)",
                 (disease_id, trait["name"]))
    return disease_id, 1


def rebuild_gene_diseases(conn: sqlite3.Connection, gene_id: int) -> None:
    """Vuelve a calcular qué enfermedades causa el gen a partir de sus variantes."""
    conn.execute("DELETE FROM gene_diseases WHERE gene_id = ?", (gene_id,))
    conn.execute(
        """INSERT OR IGNORE INTO gene_diseases (gene_id, disease_id, source_id)
           SELECT v.gene_id, vd.disease_id, vd.source_id FROM variant_diseases vd
           JOIN variants v ON v.id = vd.variant_id
           WHERE v.gene_id = ? AND vd.effect = 'causal'""", (gene_id,))
