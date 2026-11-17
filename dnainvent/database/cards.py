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
  ensayo de fase 4), o en curated.APPROVALS. «Actúan sobre la causa» si actúan sobre la
  proteína de un gen de la enfermedad, son terapia génica o celular o están en curated.CAUSE;
  si no, «otros tratamientos». Si ChEMBL registra la enfermedad con otro nombre
  (curated.DISEASE_EQUIVALENTS), se enseñan también esos fármacos, diciendo con qué nombre.

* «Retirado»: ChEMBL no dice dónde. Si curated.WITHDRAWALS no lo sabe, «en algún país».
"""

import re
import sqlite3

from ..errors import DNAinventError
from . import queries, texts
from .curated import APPROVALS, CAUSE, DISEASE_EQUIVALENTS, STRUCTURES, WITHDRAWALS
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


# Tipo de acción de ChEMBL -> clave de i18n (database.actions.*); el resto, tal cual
ACTION_KEYS = {"inhibitor": "inhibitor", "agonist": "agonist", "antagonist": "antagonist",
               "modulator": "modulator", "activator": "activator", "blocker": "blocker",
               "opener": "opener", "stabiliser": "stabiliser", "binding agent": "binding_agent",
               "positive allosteric modulator": "positive_modulator",
               "negative allosteric modulator": "negative_modulator",
               "cross-linking agent": "cross_linking", "substrate": "substrate",
               "releasing agent": "releasing_agent", "degrader": "degrader",
               "positive modulator": "positive_mod", "negative modulator": "negative_mod",
               "exogenous gene": "exogenous_gene", "partial agonist": "partial_agonist",
               "inverse agonist": "inverse_agonist", "antisense inhibitor": "antisense",
               "rnai inhibitor": "rnai", "chelating agent": "chelating", "disrupting agent": "disrupting",
               "hydrolytic enzyme": "hydrolytic_enzyme", "vaccine antigen": "vaccine_antigen",
               "exogenous protein": "exogenous_protein", "gene editing": "gene_editing"}

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


def significance_parts(value: str | None) -> list[dict]:
    """Clasificación de ClinVar en trozos traducibles. ClinVar a veces junta varias
    («pathogenic; other», «pathogenic/likely pathogenic/pathogenic, low penetrance»):
    [{'key': 'pathogenic'}, {'sep': '; '}, {'key': 'other'}]. Lo desconocido va en 'text'."""
    if not value:
        return [{"key": "none"}]
    if value in SIGNIFICANCE_KEYS:
        return [{"key": SIGNIFICANCE_KEYS[value]}]
    out = []
    for i, part in enumerate(re.split(r"(; |/)", value)):
        if i % 2:
            out.append({"sep": part})
        else:
            key = SIGNIFICANCE_KEYS.get(part)
            out.append({"key": key} if key else {"text": part})
    return out


def action_key(action: str | None) -> str | None:
    return ACTION_KEYS.get((action or "").lower())


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
    """{'en': nombre de ClinVar, 'es': ..., 'ca': ..., 'label_en': ...}. La traducción viene de la
    base de datos o de los textos escritos a mano (texts/*.json); si no hay, None y la consola
    pone el inglés con la marca «(nombre en inglés)». 'en' es siempre el de ClinVar (para
    ordenar y buscar); 'label_en' es el escrito a mano, más legible, si lo hay."""
    name = disease["name_en"]
    return {"en": name,
            "es": disease["name_es"] or texts.disease_name(name, "es"),
            "ca": disease["name_ca"] or texts.disease_name(name, "ca"),
            "label_en": texts.disease_name(name, "en")}


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
            "withdrawal": withdrawal(drug), "mechanisms": drug["mechanisms"], "cause": cause,
            "registered_as": drug["registered_as"]}


def is_approved_for(chembl_id: str | None, disease_name_en: str, max_phase, evidence) -> bool:
    """Aprobado para esa enfermedad: ficha técnica según ChEMBL, o curated.APPROVALS."""
    return (max_phase == 4 and evidence == "label") \
        or disease_name_en in APPROVALS.get(chembl_id or "", {})


def disease_drug_rows(conn: sqlite3.Connection, disease) -> list[dict]:
    """Los fármacos de la enfermedad y, si ChEMBL la registra con otro nombre, también los de
    ese nombre. Cada fármaco una vez, con su mejor fase; 'registered_as' es None o el nombre
    (disease_name) con el que está en ChEMBL."""
    rows: dict[int, dict] = {}
    names = [disease["name_en"]] + DISEASE_EQUIVALENTS.get(disease["name_en"], [])
    for position, name in enumerate(names):
        other = disease if position == 0 else conn.execute(
            "SELECT * FROM diseases WHERE name_en = ?", (name,)).fetchone()
        if other is None:
            continue
        for drug in queries.disease_drugs(conn, other["id"]):
            entry = dict(drug)
            entry["approved"] = is_approved_for(drug["chembl_id"], name, drug["phase_for_disease"],
                                                drug["evidence"])
            entry["registered_as"] = None if position == 0 else disease_name(other)
            rank = (entry["approved"], entry["phase_for_disease"] or 0, position == 0)
            current = rows.get(drug["id"])
            if current is None or rank > current["_rank"]:
                entry["_rank"] = rank
                rows[drug["id"]] = entry
    return list(rows.values())


def treatment_status(drugs: list[dict]) -> str:
    """approved · withdrawn (los aprobados se retiraron) · trials · none. Como
    queries.disease_status, pero con lo escrito a mano en curated."""
    approved = [d for d in drugs if d["approved"]]
    if any(withdrawal(d) is None for d in approved):
        return "approved"
    if approved:
        return "withdrawn"
    if any((d["phase_for_disease"] or 0) > 0 for d in drugs):
        return "trials"
    return "none"


def disease_card(conn: sqlite3.Connection, disease_id: int) -> dict:
    disease = conn.execute("SELECT * FROM diseases WHERE id = ?", (disease_id,)).fetchone()
    genes = [dict(g) for g in queries.genes_for_disease(conn, disease_id) if g["causal"]]
    gene_ids = {g["id"] for g in genes}

    drugs = disease_drug_rows(conn, disease)
    approved_cause, approved_symptoms, trials = [], [], []
    for drug in drugs:
        targets = {r[0] for r in conn.execute("SELECT gene_id FROM drug_targets WHERE drug_id = ?",
                                              (drug["id"],))}
        actions = " ".join((r[0] or "").lower() for r in conn.execute(
            "SELECT action FROM drug_mechanisms WHERE drug_id = ?", (drug["id"],)))
        cause = bool(targets & gene_ids) or (drug["molecule_type"] or "").lower() in CAUSE_MOLECULE_TYPES \
            or "gene" in actions or drug["chembl_id"] in CAUSE
        if drug["approved"]:
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
        "prevalence": prevalence(disease), "status": treatment_status(drugs),
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
    structure = STRUCTURES.get(drug["chembl_id"] or "")
    return {
        "id": drug["id"], "name": drug["name"], "chembl_id": drug["chembl_id"],
        "molecule_type": drug["molecule_type"], "first_approval": drug["first_approval"],
        "max_phase": drug["max_phase"], "withdrawal": withdrawal(drug),
        "mechanisms": unique, "genes": genes,
        "complex": any("COMPLEX" in (m["target_type"] or "") for m in unique),
        "diseases": sorted(({"id": d["id"], "name": disease_name(d),
                             "approved": is_approved_for(drug["chembl_id"], d["name_en"],
                                                         d["max_phase"], d["evidence"]),
                             "phase": d["max_phase"]} for d in diseases),
                           key=lambda d: not d["approved"]),
        "structure": structure["pdb"] if structure else None,
        "ligand": structure.get("ligand") if structure else None,
    }


# ---------------------------------------------------------------------- enlaces
# Las páginas de las fuentes de cada ficha, para ir a ver los datos originales.
NCBI = "https://www.ncbi.nlm.nih.gov"


def gene_links(card: dict) -> list[tuple[str, str]]:
    out = []
    if card.get("ncbi_gene_id"):
        out.append(("NCBI Gene", f"{NCBI}/gene/{card['ncbi_gene_id']}"))
    if card.get("uniprot_id"):
        out.append(("UniProt", f"https://www.uniprot.org/uniprotkb/{card['uniprot_id']}"))
        out.append(("AlphaFold DB", f"https://alphafold.ebi.ac.uk/entry/{card['uniprot_id']}"))
    out.append(("ClinVar", f"{NCBI}/clinvar/?term={card['symbol']}%5Bgene%5D"))
    return out


def variant_links(card: dict) -> list[tuple[str, str]]:
    out = []
    if card.get("clinvar_id"):
        out.append(("ClinVar", f"{NCBI}/clinvar/variation/{card['clinvar_id']}/"))
    if card.get("dbsnp"):
        out.append(("dbSNP", f"{NCBI}/snp/{card['dbsnp']}"))
    return out


def disease_links(card: dict) -> list[tuple[str, str]]:
    out = []
    if card.get("mondo_id"):
        out.append(("MONDO", f"https://monarchinitiative.org/{card['mondo_id']}"))
    for code in card.get("orpha_codes", [])[:2]:
        out.append((f"Orphanet {code}", f"https://www.orpha.net/en/disease/detail/{code}"))
    return out


def drug_links(card: dict) -> list[tuple[str, str]]:
    out = []
    if card.get("chembl_id"):
        out.append(("ChEMBL", f"https://www.ebi.ac.uk/chembl/compound_report_card/{card['chembl_id']}/"))
    if card.get("structure"):
        out.append(("RCSB PDB", f"https://www.rcsb.org/structure/{card['structure']}"))
    return out
