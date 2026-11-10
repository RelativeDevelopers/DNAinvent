"""Diagnóstico del recorrido de la anemia falciforme con las respuestas REALES de los servicios.

    python tools/chembl_check.py

Hace las mismas peticiones que la opción 9 (UniProt, ChEMBL y RCSB PDB), sin caché y sin tocar
tu base de datos, y enseña lo que llega: proteínas de HBB y HBA1, sus dianas en ChEMBL, los
fármacos que actúan sobre ellas, los indicados para la anemia falciforme (por MONDO y por MeSH)
y la estructura 5E83 con el voxelotor dentro. Guarda cada respuesta en la carpeta chembl_check/
(en la carpeta desde la que lo ejecutes).
"""

import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dnainvent import DNAinventError  # noqa: E402
from dnainvent.database.chembl import ChEMBLProvider, evidence, phase  # noqa: E402
from dnainvent.protein import ligands  # noqa: E402
from dnainvent.protein.cache import Cache  # noqa: E402
from dnainvent.protein.clients import http_get  # noqa: E402
from dnainvent.protein.clients.rcsb import RCSB  # noqa: E402

OUT = Path("chembl_check")
SICKLE_MONDO, SICKLE_MESH = "MONDO:0011382", "D000755"     # como lo guarda ClinVar


class Recorder:
    """Descarga como siempre y guarda cada respuesta numerada."""

    def __init__(self):
        self.n = 0
        OUT.mkdir(exist_ok=True)

    def __call__(self, url: str) -> bytes:
        raw = http_get(url, timeout=90)        # ChEMBL puede tardar mucho
        self.n += 1
        name = url.split("?")[0].rstrip("/").rsplit("/", 1)[-1].replace(".json", "")
        (OUT / f"{self.n:02d}_{name}{'' if name.endswith('.pdb') else '.json'}").write_bytes(raw)
        (OUT / "urls.txt").open("a", encoding="utf-8").write(f"{self.n:02d} {url}\n")
        return raw


def main():
    recorder = Recorder()
    (OUT / "urls.txt").write_text("", encoding="utf-8")
    with tempfile.TemporaryDirectory() as empty:
        chembl = ChEMBLProvider(http=recorder, cache=Cache(empty))
        chembl.progress = lambda key, data: print(f"   (ChEMBL va lento: intento {data['attempt']} de "
                                                  f"{data['max']} dentro de {data['wait']} s)")
        targets = {}
        print("=== Proteínas y dianas ===")
        for symbol in ("HBB", "HBA1"):
            accession = chembl.uniprot.accession(symbol)
            found = chembl.targets_for_accession(accession) if accession else []
            print(f"{symbol}: UniProt {accession} · {len(found)} dianas en ChEMBL")
            for t in found:
                accs = [c.get("accession") for c in t.get("target_components") or []]
                print(f"   {t.get('target_chembl_id')} · {t.get('pref_name')} · {t.get('target_type')} · {accs}")
                targets[t.get("target_chembl_id")] = t

        mechanisms = chembl.mechanisms(target_ids=list(targets)) if targets else []
        print(f"\n=== Mecanismos sobre esas dianas: {len(mechanisms)} ===")
        for m in mechanisms[:15]:
            print(f"   {m.get('molecule_chembl_id')} (padre {m.get('parent_molecule_chembl_id')}) · "
                  f"{m.get('action_type')} · {m.get('mechanism_of_action')} · {m.get('target_chembl_id')}")
        if mechanisms:
            print(f"   campos: {sorted(mechanisms[0])}")

        print("\n=== Indicaciones para la anemia falciforme ===")
        by_mondo = chembl.indications(efo_id=SICKLE_MONDO)
        by_mesh = chembl.indications(mesh_id=SICKLE_MESH)
        print(f"   por {SICKLE_MONDO}: {len(by_mondo)} · por MeSH {SICKLE_MESH}: {len(by_mesh)}")
        for ind in (by_mondo or by_mesh)[:15]:
            print(f"   {ind.get('parent_molecule_chembl_id') or ind.get('molecule_chembl_id')} · "
                  f"{ind.get('efo_id')} {ind.get('efo_term')!r} · {ind.get('mesh_id')} · "
                  f"fase {ind.get('max_phase_for_ind')}")
        if by_mondo or by_mesh:
            print(f"   campos: {sorted((by_mondo or by_mesh)[0])}")
            from collections import Counter
            refs = Counter(r.get("ref_type") for i in by_mondo + by_mesh for r in i.get("indication_refs") or [])
            print(f"   de dónde salen (ref_type): {dict(refs)}")
            labelled = sorted({(i.get("parent_molecule_chembl_id") or i.get("molecule_chembl_id"), i.get("max_phase_for_ind"))
                               for i in by_mondo + by_mesh if evidence(i) == "label"})
            print(f"   con ficha técnica (aprobados para ELLA si fase 4): {labelled}")

        ids = {m.get("parent_molecule_chembl_id") or m.get("molecule_chembl_id") for m in mechanisms}
        ids |= {i.get("parent_molecule_chembl_id") or i.get("molecule_chembl_id") for i in by_mondo + by_mesh}
        ids.discard(None)
        molecules = chembl.molecules(ids)
        withdrawn = chembl.withdrawn(ids)
        print(f"\n=== Fármacos: {len(molecules)} (retirados según drug_warning: {len(withdrawn)}) ===")
        approved_for_it = {i.get("parent_molecule_chembl_id") or i.get("molecule_chembl_id")
                           for i in by_mondo + by_mesh if evidence(i) == "label"}
        for m in sorted(molecules, key=lambda m: m.get("molecule_chembl_id") not in approved_for_it)[:20]:
            print(f"   {m.get('molecule_chembl_id')} · {m.get('pref_name')} · {m.get('molecule_type')} · "
                  f"fase {phase(m.get('max_phase'))} · aprobado {m.get('first_approval')} · "
                  f"withdrawn_flag={m.get('withdrawn_flag')!r} · retirado={m.get('molecule_chembl_id') in withdrawn}")
        if molecules:
            print(f"   campos: {sorted(molecules[0])}")

    print("\n=== Estructura 5E83 (hemoglobina S con voxelotor) ===")
    try:
        structure = RCSB(cache=Cache(OUT / "cache"), http=recorder).structure("5E83")
    except DNAinventError as error:
        print(f"   ERROR {error.code} {error.params}")
        return
    print(f"   {structure.title}")
    print(f"   cadenas: {structure.chain_accessions}")
    found = ligands.find_ligands(structure, include_additives=True)
    for lig in found:
        print(f"   {lig.key} · {lig.kind} · {len(lig.heavy_atoms)} átomos pesados")
    drugs = [lig for lig in found if lig.kind == "drug"]
    for drug in drugs[:2]:
        print(f"\n   Bolsillo de {drug.key}:")
        for r in ligands.pocket(structure, drug)[:10]:
            print(f"      {r.resn}{r.resi} ({r.chain}) · {r.kind} · {r.distance} Å · {r.contacts} pares")
    print(f"\nRespuestas guardadas en {OUT.resolve()}")


if __name__ == "__main__":
    try:
        main()
    except DNAinventError as error:
        print(f"\nERROR {error.code} {error.params} (las respuestas que llegaron están en {OUT.resolve()})")
