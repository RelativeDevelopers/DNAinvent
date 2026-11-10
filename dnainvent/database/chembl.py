"""ChEMBL (EMBL-EBI): fármacos, sobre qué actúan, para qué enfermedad y en qué fase están.

    provider = ChEMBLProvider()
    provider.import_gene_drugs(conn, gene_id)          # gen -> su proteína -> fármacos
    provider.import_disease_drugs(conn, disease_id)    # enfermedad -> fármacos indicados

Hacen falta los dos caminos. En la anemia falciforme:
  * el voxelotor actúa SOBRE la hemoglobina: sale desde el gen (HBA1/HBB);
  * la hidroxiurea la TRATA, pero actuando sobre otra proteína (la ribonucleótido reductasa),
    y así aumenta la hemoglobina fetal: solo sale desde la enfermedad.

API REST de ChEMBL (https://www.ebi.ac.uk/chembl/api/data/<recurso>.json), con filtros al
estilo de Django (campo__in=a,b,c) y páginas (page_meta.total_count):

    target          dianas que contienen una proteína de UniProt (proteína sola o complejo)
    mechanism       fármaco -> diana, tipo de acción (INHIBITOR...) y mecanismo en palabras
    molecule        nombre, tipo (Small molecule, Protein, Antibody, Cell...), fase máxima,
                    año de aprobación y estructura (SMILES, InChIKey)
    drug_indication fármaco -> enfermedad (EFO/MONDO y MeSH) con la fase para ESA enfermedad
    drug_warning    avisos; "Withdrawn" = retirado del mercado en algún país

Fase: 4 aprobado · 3, 2, 1 en ensayo clínico · 0.5 fase 1 temprana · sin fase: preclínico.
Que un fármaco esté aprobado para una enfermedad significa que la TRATA, no que la cure.

OJO: en drug_indication, fase 4 NO siempre es «aprobado para esta enfermedad». También sale
fase 4 cuando un fármaco ya comercializado (para otra cosa) se prueba en un ensayo de fase 4
con estos pacientes: en la anemia falciforme aparecen así estatinas, metformina u omeprazol.
Por eso cada indicación guarda de dónde viene (indication_refs): de la ficha técnica (FDA,
EMA, DailyMed, ATC) -> evidence 'label'; solo de ensayos (ClinicalTrials) -> 'trials'. Solo
'label' con fase 4 cuenta como tratamiento aprobado para esa enfermedad.

Los IDs de enfermedad de drug_indication vienen con dos puntos: efo_id = "MONDO:0011382".

Las respuestas se guardan en la caché de disco 30 días (~/.dnainvent/cache/chembl).
Datos de ChEMBL: CC BY-SA 3.0 (Zdrazil et al. 2024, Nucleic Acids Research).
"""

import json
import sqlite3
import time
import urllib.error
import urllib.parse
from dataclasses import dataclass
from typing import Callable

from ..errors import ServiceError
from ..protein.cache import Cache, key_for
from ..protein.clients import http_get
from ..protein.clients.uniprot import UniProt
from . import get_or_create_source

SERVICE = "ChEMBL"
API_URL = "https://www.ebi.ac.uk/chembl/api/data/"
PAGE_SIZE = 1000
IN_BATCH = 50                 # IDs por filtro __in (la URL no puede ser enorme)
MAX_AGE = 30 * 24 * 3600
MIN_INTERVAL = 0.2            # s entre peticiones: el EBI pide un uso razonable
# ChEMBL puede ir muy lento: comprobado en noviembre de 2026, las búsquedas por enfermedad
# tardaban ~25 s cuando respondían, y a ratos daban 500 o no respondían en 30 s. Por eso se
# espera más que en el resto de servicios y se reintenta más veces, cada vez esperando más.
TIMEOUT = 90                  # s por petición
MAX_ATTEMPTS = 5
RETRY_WAITS = (5, 10, 20, 40)  # s antes de cada reintento


def _slow_http_get(url: str) -> bytes:
    return http_get(url, timeout=TIMEOUT)

Progress = Callable[[str, dict], None]

# molecule_type de ChEMBL -> drugs.drug_type
DRUG_TYPES = {"small molecule": "small_molecule", "protein": "protein", "antibody": "antibody",
              "oligonucleotide": "oligonucleotide", "oligosaccharide": "other", "enzyme": "protein",
              "cell": "other", "gene": "other", "unknown": "other"}


def _no_progress(key: str, data: dict) -> None:
    pass


def phase(value) -> float | None:
    """'4.0', 4, 0.5... -> número; None si no hay o es -1 (desconocida)."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if number < 0 else (int(number) if number.is_integer() else number)


# Fuentes de ficha técnica. Comprobado con la anemia falciforme (nov. 2026): llegan
# ClinicalTrials, DailyMed, EMA, FDA y USAN. USAN solo registra el nombre: no es una aprobación.
LABEL_SOURCES = frozenset({"FDA", "EMA", "DailyMed", "ATC", "PMDA", "BNF", "HMA"})


def efo_to_curie(efo_id: str | None) -> str | None:
    """'MONDO_0011382' o 'MONDO:0011382' -> 'MONDO:0011382' (como los IDs de MONDO de ClinVar)."""
    return efo_id.replace("_", ":", 1) if efo_id else None


def curie_to_efo(curie: str | None) -> str | None:
    return curie.replace(":", "_", 1) if curie else None


def evidence(indication: dict) -> str | None:
    """De dónde sale una indicación: 'label' (ficha técnica de un fármaco aprobado para ella),
    'trials' (solo ensayos clínicos) o 'other'; None si ChEMBL no da referencias."""
    types = {ref.get("ref_type") for ref in indication.get("indication_refs") or []}
    types.discard(None)
    if types & LABEL_SOURCES:
        return "label"
    if "ClinicalTrials" in types:
        return "trials"
    return "other" if types else None


@dataclass
class DrugStats:
    accession: str | None = None     # proteína de UniProt del gen
    targets: int = 0                 # dianas de ChEMBL que la contienen
    drugs: int = 0                   # fármacos encontrados
    approved: int = 0                # aprobados (fase 4) para algo
    withdrawn: int = 0
    indications_linked: int = 0      # indicaciones enlazadas a enfermedades de la base de datos
    indications_other: int = 0       # indicaciones para enfermedades que no están guardadas


class ChEMBLProvider:
    service = SERVICE

    def __init__(self, http: Callable[[str], bytes] = _slow_http_get, cache: Cache | None = None,
                 sleep: Callable[[float], None] = time.sleep, uniprot: UniProt | None = None):
        self.http = http
        self.cache = cache if cache is not None else Cache()
        self.sleep = sleep
        self.uniprot = uniprot or UniProt(cache=self.cache,
                                          http=http_get if http is _slow_http_get else http)
        self._last_request = 0.0
        self.progress: Progress = _no_progress   # para avisar de los reintentos

    # ------------------------------------------------------------------ peticiones
    def get(self, resource: str, params: dict) -> dict:
        url = API_URL + resource + ".json?" + urllib.parse.urlencode(sorted(params.items()))
        key = key_for(url) + ".json"
        raw = self.cache.get("chembl", key, max_age=MAX_AGE)
        if raw is None:
            wait = MIN_INTERVAL - (time.monotonic() - self._last_request)
            if wait > 0:
                self.sleep(wait)
            self._last_request = time.monotonic()
            for attempt in range(1, MAX_ATTEMPTS + 1):
                try:
                    raw = self.http(url)
                    break
                except urllib.error.HTTPError as error:
                    error.close()
                    if error.code != 429 and error.code < 500 or attempt == MAX_ATTEMPTS:
                        raise ServiceError("service_unavailable", service=SERVICE) from None
                except (urllib.error.URLError, OSError):     # incluye el tiempo agotado
                    if attempt == MAX_ATTEMPTS:
                        raise ServiceError("service_unavailable", service=SERVICE) from None
                wait = RETRY_WAITS[attempt - 1]
                self.progress("database.stage.chembl_retry",
                              {"attempt": attempt + 1, "max": MAX_ATTEMPTS, "wait": wait})
                self.sleep(wait)
            self.cache.put("chembl", key, raw)
        try:
            data = json.loads(raw)
        except ValueError:
            raise ServiceError("service_unavailable", service=SERVICE) from None
        if not isinstance(data, dict):
            raise ServiceError("search_failed", service=SERVICE)
        return data

    def all(self, resource: str, list_key: str, params: dict) -> list[dict]:
        """Todas las páginas de un listado."""
        items: list[dict] = []
        offset = 0
        while True:
            data = self.get(resource, dict(params, limit=PAGE_SIZE, offset=offset))
            page = data.get(list_key) or []
            items.extend(page)
            total = (data.get("page_meta") or {}).get("total_count") or 0
            offset += len(page)
            if not page or offset >= total:
                return items

    def all_in(self, resource: str, list_key: str, field: str, values, params: dict | None = None) -> list[dict]:
        """Listado filtrado por una lista de IDs (campo__in), en tandas."""
        values = sorted({v for v in values if v})
        items = []
        for i in range(0, len(values), IN_BATCH):
            items += self.all(resource, list_key,
                              dict(params or {}, **{f"{field}__in": ",".join(values[i:i + IN_BATCH])}))
        return items

    # ------------------------------------------------------------------ consultas
    def targets_for_accession(self, accession: str) -> list[dict]:
        return self.all("target", "targets", {"target_components__accession": accession})

    def targets(self, target_ids) -> list[dict]:
        return self.all_in("target", "targets", "target_chembl_id", target_ids)

    def mechanisms(self, *, target_ids=None, molecule_ids=None) -> list[dict]:
        if target_ids is not None:
            return self.all_in("mechanism", "mechanisms", "target_chembl_id", target_ids)
        return self.all_in("mechanism", "mechanisms", "molecule_chembl_id", molecule_ids or [])

    def molecules(self, molecule_ids) -> list[dict]:
        return self.all_in("molecule", "molecules", "molecule_chembl_id", molecule_ids)

    def indications(self, *, molecule_ids=None, efo_id=None, mesh_id=None) -> list[dict]:
        if molecule_ids is not None:
            return self.all_in("drug_indication", "drug_indications", "molecule_chembl_id", molecule_ids)
        # Cada forma de buscar se prueba aparte: si una falla (ChEMBL da error con algunos
        # filtros) se sigue con la siguiente, y solo se avisa si fallan todas.
        attempts = []
        if efo_id:
            attempts.append({"efo_id": efo_id})
            if ":" in efo_id:                    # por si alguna versión lo escribe con "_"
                attempts.append({"efo_id": curie_to_efo(efo_id)})
        if mesh_id:
            attempts.append({"mesh_id": mesh_id})
        found, failed = [], 0
        for params in attempts:
            if found and "efo_id" in params and params["efo_id"] != efo_id:
                continue                         # el ID con ":" ya ha dado resultados
            try:
                found += self.all("drug_indication", "drug_indications", params)
            except ServiceError:
                failed += 1
        if attempts and failed == len(attempts):
            raise ServiceError("service_unavailable", service=SERVICE)
        return found

    def withdrawn(self, molecule_ids) -> set[str]:
        warnings = self.all_in("drug_warning", "drug_warnings", "molecule_chembl_id", molecule_ids)
        return {w.get("molecule_chembl_id") for w in warnings
                if (w.get("warning_type") or "").lower() == "withdrawn"}

    # ------------------------------------------------------------------ importar
    def import_gene_drugs(self, conn: sqlite3.Connection, gene_id: int,
                          progress: Progress = _no_progress) -> DrugStats:
        """Fármacos que actúan sobre la proteína del gen, y para qué enfermedades están."""
        self.progress = progress
        gene = conn.execute("SELECT * FROM genes WHERE id = ?", (gene_id,)).fetchone()
        stats = DrugStats(accession=gene["uniprot_id"])
        progress("database.stage.drug_targets", {"symbol": gene["symbol"]})
        if not stats.accession:
            stats.accession = self.uniprot.accession(gene["symbol"], gene["taxon_id"])
            if not stats.accession:
                return stats
            with conn:
                conn.execute("UPDATE genes SET uniprot_id = ? WHERE id = ?", (stats.accession, gene_id))

        targets = {t["target_chembl_id"]: t for t in self.targets_for_accession(stats.accession)
                   if t.get("target_chembl_id")}
        stats.targets = len(targets)
        mechanisms = self.mechanisms(target_ids=list(targets)) if targets else []
        self._save(conn, mechanisms, targets, progress, stats)
        return stats

    def import_disease_drugs(self, conn: sqlite3.Connection, disease_id: int,
                             progress: Progress = _no_progress) -> DrugStats:
        """Fármacos indicados para una enfermedad (aunque actúen sobre otra proteína)."""
        self.progress = progress
        disease = conn.execute("SELECT * FROM diseases WHERE id = ?", (disease_id,)).fetchone()
        stats = DrugStats()
        progress("database.stage.disease_drugs", {"disease": disease["name_en"]})
        indications = self.indications(efo_id=disease["mondo_id"], mesh_id=disease["mesh_id"])
        molecule_ids = {_molecule_id(i) for i in indications}
        if not molecule_ids:
            return stats
        mechanisms = self.mechanisms(molecule_ids=molecule_ids)
        targets = {t["target_chembl_id"]: t for t in
                   self.targets({m.get("target_chembl_id") for m in mechanisms}) if t.get("target_chembl_id")}
        stats.targets = len(targets)
        self._save(conn, mechanisms, targets, progress, stats, extra_molecules=molecule_ids,
                   indications=indications)
        return stats

    def _save(self, conn, mechanisms, targets, progress, stats, extra_molecules=(), indications=None):
        molecule_ids = {_molecule_id(m) for m in mechanisms} | set(extra_molecules)
        molecule_ids.discard(None)
        stats.drugs = len(molecule_ids)
        if not molecule_ids:
            return
        progress("database.stage.drugs", {"count": len(molecule_ids)})
        molecules = {m["molecule_chembl_id"]: m for m in self.molecules(molecule_ids)}
        withdrawn = self.withdrawn(molecule_ids)
        progress("database.stage.indications", {"count": len(molecule_ids)})
        if indications is None:
            indications = self.indications(molecule_ids=molecule_ids)
        else:   # las de la enfermedad, más las demás indicaciones de esos fármacos
            indications = indications + self.indications(molecule_ids=molecule_ids)

        with conn:
            source = get_or_create_source(conn, "chembl")
            genes_by_accession = {r["uniprot_id"]: r["id"] for r in
                                  conn.execute("SELECT id, uniprot_id FROM genes WHERE uniprot_id IS NOT NULL")}
            drug_ids = {}
            for chembl_id in molecule_ids:
                molecule = molecules.get(chembl_id, {"molecule_chembl_id": chembl_id})
                drug_ids[chembl_id] = _upsert_drug(conn, molecule, chembl_id in withdrawn, source)
                stats.approved += phase(molecule.get("max_phase")) == 4
                stats.withdrawn += chembl_id in withdrawn

            for m in mechanisms:
                drug_id = drug_ids.get(_molecule_id(m))
                if drug_id is None:
                    continue
                target = targets.get(m.get("target_chembl_id"), {})
                action = (m.get("action_type") or "").lower() or None
                mechanism = m.get("mechanism_of_action") or ""
                conn.execute(
                    """INSERT OR REPLACE INTO drug_mechanisms
                       (drug_id, target_chembl_id, target_name, target_type, action, mechanism)
                       VALUES (?, ?, ?, ?, ?, ?)""",
                    (drug_id, m.get("target_chembl_id") or "", target.get("pref_name"),
                     target.get("target_type"), action, mechanism))
                for component in target.get("target_components") or []:
                    gene_id = genes_by_accession.get(component.get("accession"))
                    if gene_id is not None:
                        conn.execute(
                            """INSERT INTO drug_targets (drug_id, gene_id, action, mechanism, source_id)
                               VALUES (?, ?, ?, ?, ?)
                               ON CONFLICT (drug_id, gene_id) DO UPDATE SET action = excluded.action,
                               mechanism = excluded.mechanism""",
                            (drug_id, gene_id, action, mechanism, source))

            for ind in indications:
                drug_id = drug_ids.get(_molecule_id(ind))
                if drug_id is None:
                    continue
                disease_id = _match_disease(conn, ind)
                if disease_id is None:
                    stats.indications_other += 1
                    continue
                stats.indications_linked += 1
                # La misma pareja fármaco-enfermedad puede llegar dos veces (por MONDO y por MeSH):
                # se queda la fase más alta y, si alguna viene de la ficha técnica, 'label'
                conn.execute(
                    """INSERT INTO drug_diseases (drug_id, disease_id, max_phase, source_id, indication,
                                                  efo_id, evidence)
                       VALUES (?, ?, ?, ?, ?, ?, ?)
                       ON CONFLICT (drug_id, disease_id) DO UPDATE SET
                         max_phase = MAX(COALESCE(drug_diseases.max_phase, -1), COALESCE(excluded.max_phase, -1)),
                         indication = excluded.indication, efo_id = excluded.efo_id,
                         evidence = CASE WHEN 'label' IN (drug_diseases.evidence, excluded.evidence)
                                         THEN 'label' ELSE COALESCE(excluded.evidence, drug_diseases.evidence) END""",
                    (drug_id, disease_id, phase(ind.get("max_phase_for_ind")), source,
                     ind.get("efo_term") or ind.get("mesh_heading"), ind.get("efo_id"), evidence(ind)))


# ---------------------------------------------------------------------- escritura
def _molecule_id(item: dict) -> str | None:
    """El fármaco «padre»: la sal (imatinib mesilato) y la base (imatinib) son el mismo."""
    return item.get("parent_molecule_chembl_id") or item.get("molecule_chembl_id")


def _upsert_drug(conn, molecule: dict, withdrawn: bool, source_id: int) -> int:
    chembl_id = molecule["molecule_chembl_id"]
    structures = molecule.get("molecule_structures") or {}
    molecule_type = molecule.get("molecule_type")
    first_approval = molecule.get("first_approval")
    values = (molecule.get("pref_name") or chembl_id,
              DRUG_TYPES.get((molecule_type or "").lower(), "other"),
              structures.get("canonical_smiles"), structures.get("standard_inchi_key"),
              phase(molecule.get("max_phase")), molecule_type,
              int(first_approval) if first_approval else None,
              1 if withdrawn or molecule.get("withdrawn_flag") else 0, source_id)
    row = conn.execute("SELECT id FROM drugs WHERE chembl_id = ?", (chembl_id,)).fetchone()
    if row:
        conn.execute(
            """UPDATE drugs SET name = ?, drug_type = ?, smiles = ?, inchikey = ?, max_phase = ?,
               molecule_type = ?, first_approval = ?, withdrawn = ?, source_id = ? WHERE id = ?""",
            values + (row["id"],))
        return row["id"]
    drug_id = conn.execute(
        """INSERT INTO drugs (name, drug_type, smiles, inchikey, max_phase, molecule_type,
                              first_approval, withdrawn, source_id, chembl_id)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""", values + (chembl_id,)).lastrowid
    conn.execute("INSERT INTO search_index (kind, ref_id, text) VALUES ('drug', ?, ?)",
                 (drug_id, f"{values[0]} {chembl_id}"))
    return drug_id


def _match_disease(conn, indication: dict) -> int | None:
    """La enfermedad guardada que corresponde a una indicación: por MONDO, por MeSH o por nombre."""
    curie = efo_to_curie(indication.get("efo_id"))
    if curie and curie.startswith("MONDO:"):
        row = conn.execute("SELECT id FROM diseases WHERE mondo_id = ?", (curie,)).fetchone()
        if row:
            return row["id"]
    if indication.get("mesh_id"):
        row = conn.execute("SELECT id FROM diseases WHERE mesh_id = ?", (indication["mesh_id"],)).fetchone()
        if row:
            return row["id"]
    for name in (indication.get("efo_term"), indication.get("mesh_heading")):
        if name:
            row = conn.execute("SELECT id FROM diseases WHERE name_en = ? COLLATE NOCASE", (name,)).fetchone()
            if row:
                return row["id"]
    return None
