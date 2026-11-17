"""Las fichas de la opción 9: gen, mutación, enfermedad y fármaco.

    card = gene_card(conn, gene_id)
    card["colors"]            # {'pathogenic': 409, 'uncertain': 384, 'benign': 844, ...}

Devuelven datos (diccionarios y listas), sin textos traducidos ni colores de terminal: la
consola los pinta con colorama y la web (fase 2) los convertirá en HTML. Los textos que
dependen del idioma se devuelven como CLAVES (p. ej. 'autosomal_recessive') para i18n.

Reglas para no contar cosas falsas (decididas con Alexei el 17-11-2026):

* Colores de una variante, según su clasificación de ClinVar: rojo (causa enfermedad),
  amarillo (incierta o con conflictos), verde (sana), especial (protectora, factor de riesgo,
  asociación, respuesta a un fármaco) y gris (sin clasificar: «other», «not provided»...).

* Nivel de revisión: el estado de revisión de ClinVar en porcentaje. NO es la probabilidad de
  que la clasificación sea cierta, sino cuánto la han revisado.

* Enfermedades de una variante: algunos laboratorios mandan a ClinVar UN registro (RCV) con
  varias enfermedades a la vez y sin criterios; tal cual, diría que HbS causa malaria. En la
  lista principal solo cuentan los registros de UNA enfermedad con al menos un 25 % de revisión
  («fiables»); los demás se cuentan aparte, para los detalles técnicos.

* Fármacos de una enfermedad: aprobado = con ficha técnica para ESA enfermedad (no basta un
  ensayo de fase 4). «Actúan sobre la causa» si actúan sobre la proteína de un gen de la
  enfermedad o son terapia génica o celular; si no, «alivian síntomas o complicaciones».

* «Retirado»: ChEMBL no dice dónde. Si curated.WITHDRAWALS no lo sabe, «en algún país».
"""

import re
import sqlite3

from ..errors import DNAinventError
from . import queries
from .curated import STRUCTURES, WITHDRAWALS
from .variants import locate, side_by_side

COLORS = ("pathogenic", "uncertain", "benign", "special", "unclassified")

REVIEW_LEVELS = {                      # estado de revisión de ClinVar -> (porcentaje, clave)
    "practice guideline": (100, "guideline"),
    "reviewed by expert panel": (75, "expert_panel"),
    "criteria provided, multiple submitters, no conflicts": (50, "multiple"),
    "criteria provided, conflicting classifications": (25, "conflicting"),
    "criteria provided, conflicting interpretations": (25, "conflicting"),
    "criteria provided, single submitter": (25, "single"),
}
TRUSTED_REVIEW = 25                    # % mínimo para la lista principal

INHERITANCE_KEYS = {
    "autosomal recessive": "autosomal_recessive", "autosomal dominant": "autosomal_dominant",
    "x-linked recessive": "x_recessive", "x-linked dominant": "x_dominant",
    "y-linked": "y_linked", "mitochondrial inheritance": "mitochondrial",
    "multigenic/multifactorial": "multifactorial", "semi-dominant": "semi_dominant",
}
PREVALENCE_TYPE_KEYS = {"point prevalence": "point", "birth prevalence": "birth",
                        "annual incidence": "annual", "lifetime prevalence": "lifetime",
                        "lifelong prevalence": "lifetime"}
PLACE_KEYS = {"worldwide": "worldwide", "europe": "europe"}

# Nombres de ClinVar que no son una enfermedad concreta: van al final de las listas
GENERIC_DISEASE = re.compile(r"(-related (disorder|condition)s?$|^hereditary cancer-predisposing "
                             r"syndrome$|^cardiovascular phenotype$|^inborn genetic diseases$|"
                             r"^see cases$|cancer predisposition$|diagnostic test$)", re.IGNORECASE)

CAUSE_MOLECULE_TYPES = {"gene", "cell"}


# ---------------------------------------------------------------------- piezas
def color(significance: str | None) -> str:
    """Clasificación de ClinVar -> uno de COLORS. Si junta varias, manda la más grave
    («pathogenic; other» es roja), salvo los conflictos, que siempre son inciertos."""
    s = (significance or "").lower()
    if "conflicting" in s:
        return "uncertain"
    if "pathogenic" in s and "benign" not in s:
        return "pathogenic"
    if "uncertain" in s:
        return "uncertain"
    if "benign" in s and "pathogenic" not in s:
        return "benign"
    if any(word in s for word in ("protective", "risk factor", "risk allele", "association",
                                  "drug response", "affects")):
        return "special"
    return "unclassified"


def review_level(status: str | None) -> tuple[int, str]:
    """(porcentaje, clave). Lo desconocido o «sin criterios» es 0 %."""
    s = (status or "").strip().lower()
    if s in REVIEW_LEVELS:
        return REVIEW_LEVELS[s]
    return (0, "no_criteria" if s.startswith("no ") or not s else "other")


def frequency(value: float | None) -> dict | None:
    """Frecuencia alélica -> {'value': 0.01298, 'one_in': 77} (1 de cada 77 copias del gen)."""
    if value is None or value <= 0:
        return None
    return {"value": value, "one_in": max(1, round(1 / value))}


def is_generic(name: str | None) -> bool:
    return bool(GENERIC_DISEASE.search(name or ""))


def inheritance(text: str | None) -> list[str]:
    """'Autosomal recessive|X-linked recessive' -> ['autosomal_recessive', 'x_recessive']
    (los que no se conocen se devuelven tal cual)."""
    return [INHERITANCE_KEYS.get(p.strip().lower(), p.strip()) for p in (text or "").split("|") if p.strip()]


def prevalence(disease) -> dict | None:
    """{'class': '1-5 / 10 000', 'type': 'point', 'place': 'europe'} o None."""
    if not disease["prevalence_class"]:
        return None
    place = (disease["prevalence_geo"] or "").lower()
    return {"class": disease["prevalence_class"],
            "type": PREVALENCE_TYPE_KEYS.get((disease["prevalence_type"] or "").lower(), "other"),
            "place": PLACE_KEYS.get(place, disease["prevalence_geo"] or "")}


def withdrawal(drug) -> dict | None:
    """None si no está retirado; si no, {'scope': 'worldwide'|'some', 'year': ..., 'where': ...}."""
    known = WITHDRAWALS.get(drug["chembl_id"] or "")
    if known:
        return dict(known)
    if drug["withdrawn"]:
        return {"scope": "some"}
    return None


def disease_name(disease) -> dict:
    """{'en': ..., 'es': ..., 'ca': ...}: los que no estén traducidos son None (la consola pone
    el inglés con la marca «(nombre en inglés)»)."""
    return {"en": disease["name_en"], "es": disease["name_es"], "ca": disease["name_ca"]}


def _trusted_rcvs(conn: sqlite3.Connection) -> str:
    """Subconsulta SQL: los RCV de UNA sola enfermedad."""
    return ("SELECT rcv_accession FROM variant_diseases WHERE rcv_accession IS NOT NULL "
            "GROUP BY rcv_accession HAVING COUNT(DISTINCT disease_id) = 1")


def _is_trusted(row, single: set[str]) -> bool:
    return row["rcv_accession"] in single and review_level(row["review_status"])[0] >= TRUSTED_REVIEW


# ---------------------------------------------------------------------- gen
def main_diseases(conn: sqlite3.Connection, gene_id: int, limit: int = 5) -> list[dict]:
    """Enfermedades del gen por número de variantes que la CAUSAN según registros fiables.
    Los nombres genéricos («HBB-related disorder») van detrás."""
    rows = conn.execute(
        f"""SELECT d.*, COUNT(DISTINCT vd.variant_id) AS causal
            FROM variant_diseases vd JOIN variants v ON v.id = vd.variant_id
            JOIN diseases d ON d.id = vd.disease_id
            WHERE v.gene_id = ? AND vd.effect = 'causal'
              AND vd.rcv_accession IN ({_trusted_rcvs(conn)})
              AND vd.review_status IN ({",".join("?" * len(REVIEW_LEVELS))})
            GROUP BY d.id""", (gene_id, *REVIEW_LEVELS)).fetchall()
    rows = sorted(rows, key=lambda r: (is_generic(r["name_en"]), -r["causal"], r["name_en"]))
    return [{"id": r["id"], "name": disease_name(r), "causal": r["causal"],
             "inheritance": inheritance(r["inheritance"]), "prevalence": prevalence(r),
             "generic": is_generic(r["name_en"])} for r in rows[:limit]]


def gene_card(conn: sqlite3.Connection, gene_id: int) -> dict:
    gene = conn.execute("SELECT * FROM genes WHERE id = ?", (gene_id,)).fetchone()
    counts = dict.fromkeys(COLORS, 0)
    for (significance, n) in conn.execute(
            "SELECT clinical_significance, COUNT(*) FROM variants WHERE gene_id = ? "
            "GROUP BY clinical_significance", (gene_id,)):
        counts[color(significance)] += n
    named = conn.execute("SELECT id, common_name, hgvs_c, hgvs_p FROM variants WHERE gene_id = ? "
                         "AND common_name IS NOT NULL ORDER BY id", (gene_id,)).fetchall()
    mrna = queries.reference_mrna(conn, gene_id)
    return {
        "id": gene["id"], "symbol": gene["symbol"], "name": gene["name"],
        "chromosome": gene["chromosome"], "ncbi_gene_id": gene["ncbi_gene_id"],
        "uniprot_id": gene["uniprot_id"],
        "reference_mrna": f"{mrna['accession']}.{mrna['version']}" if mrna else None,
        "colors": counts, "total": sum(counts.values()),
        "diseases": main_diseases(conn, gene_id),
        "known_variants": [dict(r) for r in named],
        "drugs": len(queries.gene_drugs(conn, gene_id)),
    }


# ---------------------------------------------------------------------- mutación
def featured_variants(conn: sqlite3.Connection, gene_id: int, limit: int = 20) -> list[sqlite3.Row]:
    """Las variantes que se enseñan sin buscar nada: primero las de nombre propio (HbS,
    ΔF508), luego las que causan enfermedad mejor revisadas y, a igualdad, las más frecuentes.
    Así salen las conocidas y no grandes deleciones sin revisar."""
    rows = conn.execute(
        """SELECT v.*, s.accession, s.version FROM variants v JOIN sequences s ON s.id = v.sequence_id
           WHERE v.gene_id = ?""", (gene_id,)).fetchall()
    chosen = [r for r in rows if r["common_name"] or color(r["clinical_significance"]) == "pathogenic"]
    chosen.sort(key=lambda r: (r["common_name"] is None, -review_level(r["review_status"])[0],
                               r["variant_type"] != "snv", -(r["population_freq"] or 0), r["hgvs_c"] or ""))
    return chosen[:limit]



def variant_card(conn: sqlite3.Connection, variant_id: int) -> dict:
    v = conn.execute(
        """SELECT v.*, g.symbol, s.accession, s.version, s.residues IS NOT NULL AS has_sequence
           FROM variants v JOIN genes g ON g.id = v.gene_id JOIN sequences s ON s.id = v.sequence_id
           WHERE v.id = ?""", (variant_id,)).fetchone()
    percent, review_key = review_level(v["review_status"])

    change = None
    if v["has_sequence"]:
        try:
            effect, residues, cds_start = locate(conn, variant_id)
        except DNAinventError:
            effect = None               # no es una sustitución en el CDS, o no coincide
        if effect is not None:
            change = {"effect": effect, "view": side_by_side(effect, residues, cds_start)}

    single = {r[0] for r in conn.execute(_trusted_rcvs(conn))}
    rows = conn.execute(
        """SELECT vd.*, d.name_en, d.name_es, d.name_ca FROM variant_diseases vd
           JOIN diseases d ON d.id = vd.disease_id WHERE vd.variant_id = ?""", (variant_id,)).fetchall()
    best: dict[int, dict] = {}
    others: set[int] = set()
    for row in rows:
        if not _is_trusted(row, single):
            others.add(row["disease_id"])
            continue
        level = review_level(row["review_status"])[0]
        current = best.get(row["disease_id"])
        if current is None or level > current["review"]:
            best[row["disease_id"]] = {
                "id": row["disease_id"], "name": disease_name(row),
                "effect": row["effect"], "color": color(row["clinical_significance"]),
                "significance": row["clinical_significance"], "review": level,
                "rcv": row["rcv_accession"]}
    diseases = sorted(best.values(), key=lambda d: (-d["review"], is_generic(d["name"]["en"]),
                                                    d["name"]["en"]))
    return {
        "id": v["id"], "gene": v["symbol"], "gene_id": v["gene_id"],
        "common_name": v["common_name"], "hgvs_c": v["hgvs_c"], "hgvs_p": v["hgvs_p"],
        "accession": f"{v['accession']}.{v['version']}", "dbsnp": v["dbsnp_rs"],
        "clinvar_id": v["clinvar_id"], "variant_type": v["variant_type"],
        "color": color(v["clinical_significance"]), "significance": v["clinical_significance"],
        "review": percent, "review_key": review_key,
        "frequency": frequency(v["population_freq"]), "frequency_source": v["freq_source"],
        "change": change,
        "diseases": diseases,
        "other_diseases": len(others - set(best)),
    }


# ---------------------------------------------------------------------- enfermedad
def _drug_entry(drug, cause: bool) -> dict:
    return {"id": drug["id"], "name": drug["name"], "chembl_id": drug["chembl_id"],
            "first_approval": drug["first_approval"], "molecule_type": drug["molecule_type"],
            "withdrawal": withdrawal(drug), "mechanisms": drug["mechanisms"], "cause": cause}


def disease_card(conn: sqlite3.Connection, disease_id: int) -> dict:
    disease = conn.execute("SELECT * FROM diseases WHERE id = ?", (disease_id,)).fetchone()
    genes = [dict(g) for g in queries.genes_for_disease(conn, disease_id) if g["causal"]]
    gene_ids = {g["id"] for g in genes}

    approved_cause, approved_symptoms, trials = [], [], []
    for drug in queries.disease_drugs(conn, disease_id):
        targets = {r[0] for r in conn.execute("SELECT gene_id FROM drug_targets WHERE drug_id = ?",
                                              (drug["id"],))}
        actions = " ".join((r[0] or "").lower() for r in conn.execute(
            "SELECT action FROM drug_mechanisms WHERE drug_id = ?", (drug["id"],)))
        cause = bool(targets & gene_ids) or (drug["molecule_type"] or "").lower() in CAUSE_MOLECULE_TYPES \
            or "gene" in actions
        if drug["phase_for_disease"] == 4 and drug["evidence"] == "label":
            (approved_cause if cause else approved_symptoms).append(_drug_entry(drug, cause))
        elif drug["phase_for_disease"] is not None and drug["phase_for_disease"] > 0:
            trials.append(_drug_entry(drug, cause))

    def order(items):          # los que siguen en el mercado primero; luego por año
        return sorted(items, key=lambda d: (d["withdrawal"] is not None, d["first_approval"] or 9999,
                                            d["name"]))

    return {
        "id": disease["id"], "name": disease_name(disease), "mondo_id": disease["mondo_id"],
        "orpha_codes": [c for c in (disease["orpha_codes"] or "").split("|") if c],
        "genes": genes, "inheritance": inheritance(disease["inheritance"]),
        "onset": [p for p in (disease["onset"] or "").split("|") if p],
        "prevalence": prevalence(disease), "status": queries.disease_status(conn, disease_id),
        "approved_cause": order(approved_cause), "approved_symptoms": order(approved_symptoms),
        "trials": sorted(trials, key=lambda d: d["name"]),
    }


# ---------------------------------------------------------------------- fármaco
def drug_card(conn: sqlite3.Connection, drug_id: int) -> dict:
    drug = conn.execute("SELECT * FROM drugs WHERE id = ?", (drug_id,)).fetchone()
    mechanisms = [dict(r) for r in conn.execute(
        "SELECT target_name, target_type, action, mechanism FROM drug_mechanisms WHERE drug_id = ? "
        "AND (target_name IS NOT NULL OR action IS NOT NULL) ORDER BY target_name", (drug_id,))]
    # Sin repetir la misma diana con la misma acción (ivacaftor: activador y modulador positivo)
    seen, unique = set(), []
    for m in mechanisms:
        key = (m["target_name"], m["action"])
        if key not in seen:
            seen.add(key)
            unique.append(m)
    genes = [dict(r) for r in queries.drug_genes(conn, drug_id)]
    diseases = conn.execute(
        """SELECT d.id, d.name_en, d.name_es, d.name_ca, dd.max_phase, dd.evidence
           FROM drug_diseases dd JOIN diseases d ON d.id = dd.disease_id WHERE dd.drug_id = ?
           ORDER BY dd.max_phase = 4 AND dd.evidence = 'label' DESC, dd.max_phase DESC, d.name_en""",
        (drug_id,)).fetchall()
    return {
        "id": drug["id"], "name": drug["name"], "chembl_id": drug["chembl_id"],
        "molecule_type": drug["molecule_type"], "first_approval": drug["first_approval"],
        "max_phase": drug["max_phase"], "withdrawal": withdrawal(drug),
        "mechanisms": unique, "genes": genes,
        "diseases": [{"id": d["id"], "name": disease_name(d),
                      "approved": d["max_phase"] == 4 and d["evidence"] == "label",
                      "phase": d["max_phase"]} for d in diseases],
        "structure": STRUCTURES.get(drug["chembl_id"] or ""),
    }
