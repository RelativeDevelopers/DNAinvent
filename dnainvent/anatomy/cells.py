"""Qué células hay en cada parte del cuerpo y qué genes las marcan.

Fuente principal: las tablas ASCT+B («Anatomical Structures, Cell Types and Biomarkers») del
Human Reference Atlas (HuBMAP), CC BY 4.0. Cada fila es una cadena de partes del cuerpo (de la
más grande a la más pequeña), los tipos de célula que hay en la última y los genes que los
marcan:

    heart > right atrium > sinoatrial node   →   cardiac pacemaker cell   →   HCN4, GJC1

Las partes del cuerpo usan identificadores de Uberon (UBERON:0000948) o de la FMA (FMA:3823);
las células, de la Cell Ontology (CL:1000477); los genes, de HGNC. BodyParts3D usa la FMA, así
que para enlazar se usa también Uberon, que dice a qué concepto de la FMA equivale cada uno
(«xref: FMA:...» en uberon-edit.obo).

Si una pieza no tiene células propias en las tablas, se sube por el árbol «forma parte de»
(un trozo de la pared del ventrículo → el corazón). Y si tampoco, se usa el tejido general de su
sistema (tissues.py): el músculo esquelético, el hueso, las arterias...
"""

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path

ASCTB_URL = ("https://raw.githubusercontent.com/hubmapconsortium/ccf-releases/main/v2.0/asct-b/"
             "asct-b-vh-{organ}.csv")
# Las tablas de la versión 2.0 del Human Reference Atlas (diciembre de 2023) que se pueden bajar
# de GitHub. Las demás (cerebro, ojo, intestinos...) no están en esta versión.
ASCTB_ORGANS = ("heart", "lung", "kidney", "liver", "skin", "bone-marrow", "blood-vasculature",
                "lymph-node", "thymus", "pancreas", "prostate", "knee", "peripheral-nervous-system",
                "skeleton", "muscular-system", "main-bronchus", "lymph-vasculature",
                "anatomical-systems")
ASCTB_CITATION = ("Human Reference Atlas, ASCT+B tables v2.0 (HuBMAP), CC BY 4.0 — "
                  "Börner et al. 2021, Nature Cell Biology")
UBERON_URL = "https://raw.githubusercontent.com/obophenotype/uberon/master/src/ontology/uberon-edit.obo"


@dataclass
class AsctbRow:
    organ: str                                   # la tabla de la que sale ('heart')
    structures: list                             # [(id, etiqueta)], de la más grande a la más pequeña
    cells: list = field(default_factory=list)    # [(CL:..., etiqueta)]
    genes: list = field(default_factory=list)    # [(símbolo, nombre, HGNC:...)]


def normalize_id(value: str) -> str:
    """'FMA:3823' y 'FMA3823' son lo mismo; BodyParts3D usa la segunda forma."""
    value = value.strip()
    return "FMA" + value[4:] if value.startswith("FMA:") else value


def read_asctb(path, organ: str | None = None) -> list[AsctbRow]:
    """Una tabla ASCT+B en CSV (las primeras líneas son autores y versión; los datos empiezan
    en la fila cuya primera columna es «AS/1»)."""
    path = Path(path)
    organ = organ or path.stem.replace("asct-b-vh-", "").replace("asctb-", "")
    with open(path, encoding="utf-8", newline="") as f:
        rows = list(csv.reader(f))
    start = next((i for i, r in enumerate(rows) if r and r[0].strip() == "AS/1"), None)
    if start is None:
        return []
    header = rows[start]

    def columns(kind: str):
        out = []
        for i, name in enumerate(header):
            m = re.fullmatch(rf"{kind}/(\d+)", name.strip())
            if m:
                n = m.group(1)
                out.append((i, header.index(f"{kind}/{n}/LABEL") if f"{kind}/{n}/LABEL" in header else None,
                            header.index(f"{kind}/{n}/ID") if f"{kind}/{n}/ID" in header else None))
        return out

    as_cols, ct_cols, gene_cols = columns("AS"), columns("CT"), columns("BGene")
    out = []
    for raw in rows[start + 1:]:
        raw = raw + [""] * (len(header) - len(raw))

        def pick(cols, ids_needed=True):
            items = []
            for name_i, label_i, id_i in cols:
                ident = raw[id_i].strip() if id_i is not None else ""
                label = (raw[label_i].strip() if label_i is not None else "") or raw[name_i].strip()
                if ident or (label and not ids_needed):
                    items.append((normalize_id(ident), label, raw[name_i].strip()))
            return items

        structures = [(i, label) for i, label, _ in pick(as_cols) if i]
        if not structures:
            continue
        cells = [(i, label) for i, label, _ in pick(ct_cols) if i.startswith("CL:")]
        genes = [(name or label, label, i) for i, label, name in pick(gene_cols) if i.startswith("HGNC:")]
        out.append(AsctbRow(organ, structures, cells, genes))
    return out


def read_uberon_fma(path) -> dict[str, set[str]]:
    """FMA -> {UBERON:...}, de las líneas «xref: FMA:...» de uberon-edit.obo."""
    out: dict[str, set[str]] = {}
    current = None
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("["):
                current = None
            elif line.startswith("id: UBERON:"):
                current = line[4:].strip()
            elif current and line.startswith("xref: FMA:"):
                fma = "FMA" + line[len("xref: FMA:"):].split()[0].strip()
                out.setdefault(fma, set()).add(current)
    return out


def _plain(text: str) -> str:
    return " ".join(text.lower().split())


class CellIndex:
    """Las células (y sus genes) de cada parte del cuerpo de las tablas ASCT+B.

    Cada fila es una cadena «corazón > aurícula derecha > nodo sinusal → célula marcapasos».
    La célula es, sobre todo, de la parte más pequeña (el nodo sinusal: «leaf»). Pero las partes
    grandes también la contienen: el corazón entero junta las células de todas sus filas
    («parts»). Al buscar se prefiere la parte exacta; las entradas sin ningún gen no cuentan."""

    def __init__(self, rows: list[AsctbRow], fma_to_uberon: dict[str, set[str]]):
        self.fma_to_uberon = fma_to_uberon
        self.by_structure: dict[str, dict] = {}       # la parte más pequeña de cada fila
        self.by_part: dict[str, dict] = {}            # cualquier parte de la fila (lo que contiene)
        for row in rows:
            if not row.cells:
                continue
            leaf_id, leaf_label = row.structures[-1]
            targets = [(self.by_structure, leaf_id, leaf_label, "leaf")]
            targets += [(self.by_part, sid, slabel, "parts") for sid, slabel in row.structures]
            for table, sid, slabel, kind in targets:
                entry = table.setdefault(sid, {"key": sid, "label": slabel, "organ": row.organ,
                                               "kind": kind, "cells": {}})
                for cl, label in row.cells:
                    cell = entry["cells"].setdefault(cl, {"label": label, "genes": {}})
                    for symbol, name, hgnc in row.genes:
                        cell["genes"][symbol] = (name, hgnc)
        # por el nombre, para cuando Uberon no da la equivalencia (el «kidney» de la FMA, 7203,
        # es «adult mammalian kidney» en Uberon, y la tabla del riñón usa «kidney»)
        for key, entry in self.by_part.items():
            leaf = self.by_structure.get(key)
            if leaf is not None and set(leaf["cells"]) == set(entry["cells"]):
                entry["kind"] = "leaf"
        self.by_name = {"parts": {_plain(e["label"]): e for e in self.by_part.values() if e["label"]}}

    @staticmethod
    def _has_genes(entry: dict) -> bool:
        return any(cell["genes"] for cell in entry["cells"].values())

    def _ids(self, fma: str) -> list[str]:
        return [fma] + sorted(self.fma_to_uberon.get(fma, ()))

    def direct(self, fma: str, name: str = "") -> tuple[str, dict] | None:
        """(identificador con el que se encontró, entrada) para un concepto de la FMA, por su ID
        (o el de Uberon) y, si no, por su nombre en inglés. La entrada junta las células de esa
        parte y las de sus partes («kind» = «parts» si las de sus partes añaden alguna). Solo
        entradas con algún gen."""
        for table_name, table in (("parts", self.by_part),):
            for key in self._ids(fma):
                entry = table.get(key)
                if entry and self._has_genes(entry):
                    return key, entry
            entry = self.by_name[table_name].get(_plain(name)) if name else None
            if entry and self._has_genes(entry):
                return entry["key"], entry
        return None

    def find(self, bp, fma: str, max_depth: int = 6):
        """Sube por «forma parte de» (y por «es un»: el pulmón derecho es un pulmón) hasta
        encontrar células: (concepto FMA donde se encontraron, cuántos pasos se subió, entrada)
        o None."""
        frontier, seen = [fma], {fma}
        for depth in range(max_depth + 1):
            for concept in sorted(frontier):
                hit = self.direct(concept, bp.name(concept))
                if hit:
                    return concept, depth, hit[1]
            nxt = []
            for concept in frontier:
                for parent in sorted(bp.partof_parents(concept)) + sorted(bp.isa_parents(concept)):
                    if parent not in seen:
                        seen.add(parent)
                        nxt.append(parent)
            if not nxt:
                return None
            frontier = nxt
        return None
