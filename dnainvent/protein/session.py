"""Sesión del microscopio para proteínas: lo que la consola le pasa al microscopio.

El microscopio se abre en otro proceso (la consola sigue funcionando), así que no se le
pueden pasar objetos de Python: se guardan las estructuras en PDB y un JSON que las describe.

    {
      "version": 1,
      "panels": [                                     una (sola) o dos (comparación)
        {"pdb": ".../AF-P68871-F1-v6.pdb", "label": "Hemoglobin subunit beta · P68871",
         "highlight": {"6": "Glu7Val"}, "start": 0, "count": 147},
        {"pdb": ".../ESMFold_5b53af62.pdb", "label": "Tu proteína (ESMFold)",
         "highlight": {"6": "Glu7Val"}, "start": 0, "count": 147}
      ],
      "origin": [x, y, z]                             Å: centro común de las dos
    }

La segunda estructura ya va superpuesta sobre la primera (compare.moved), y las dos se
centran en el mismo origen: así comparten coordenadas y la nave las ve en el mismo sitio.

Python puro (sin NumPy): la consola lo usa aunque el microscopio no esté instalado.
"""

import json
from dataclasses import dataclass
from pathlib import Path

from ..scene.limits import RECOMMENDED_RESIDUES
from .cache import Cache, key_for
from .compare import moved
from .pdb import Structure, read_pdb, write_pdb

VERSION = 1


@dataclass
class Panel:
    structure: Structure
    label: str
    highlight: dict[int, str]          # aminoácido (desde 0) -> "Glu7Val"
    start: int = 0
    count: int | None = None


def panels_for(report, labels: dict[str, str], start: int = 0,
               count: int | None = RECOMMENDED_RESIDUES) -> list[Panel]:
    """Qué enseñar de un ProteinReport (pipeline.py):

        conocida            la natural
        variante / pariente la natural y la tuya (si se pudieron obtener las dos)
        inventada           la tuya (predicha)

    labels: {"natural": "...", "yours": "..."} ya traducidos.
    start, count: tramo de la natural (o de la tuya si no hay natural); en la otra se toma
    el tramo equivalente según el alineamiento.
    """
    comparison = report.comparison
    natural_marks, yours_marks = {}, {}
    if comparison is not None:
        for d in comparison.differences:
            natural_marks[d.natural_position - 1] = d.label
            yours_marks[d.user_position - 1] = d.label

    panels = []
    if report.structure is not None:
        panels.append(Panel(report.structure, labels["natural"], natural_marks, start, count))
    if report.prediction is not None:
        prediction = report.prediction
        if comparison is not None and comparison.superposition is not None:
            prediction = moved(prediction, comparison.superposition)
        panel = Panel(prediction, labels["yours"], yours_marks, start, count)
        if panels and comparison is not None:
            panel.start, panel.count = _mapped_window(comparison.pairs, start, count,
                                                      len(prediction.residues()))
        panels.append(panel)
    return panels


def _mapped_window(pairs, start, count, length) -> tuple[int, int]:
    """Tramo de la tuya que corresponde al tramo [start, start + count) de la natural."""
    end = start + count if count is not None else float("inf")
    mine = [i for i, j in pairs if start <= j < end]
    if not mine:
        return 0, length
    return min(mine), max(mine) - min(mine) + 1


def write_session(panels: list[Panel], cache: Cache | None = None) -> Path:
    """Guarda las estructuras y el JSON en la caché (~/.dnainvent/cache/microscope)."""
    cache = cache if cache is not None else Cache()
    entries = []
    for panel in panels:
        text = write_pdb(panel.structure.atoms, panel.structure.title)
        path = cache.put("microscope", key_for(text)[:16] + ".pdb", text.encode("utf-8"))
        entries.append({"pdb": str(path), "label": panel.label,
                        "highlight": {str(k): v for k, v in panel.highlight.items()},
                        "start": panel.start, "count": panel.count})
    data = {"version": VERSION, "panels": entries, "origin": common_origin(panels)}
    text = json.dumps(data, ensure_ascii=False, indent=1)
    return cache.put("microscope", "session-" + key_for(text)[:16] + ".json",
                     text.encode("utf-8"))


def common_origin(panels: list[Panel]) -> list[float]:
    """Centro (Å) del tramo de la primera estructura: las dos se centran ahí."""
    first = panels[0]
    residues = first.structure.residues()
    end = len(residues) if first.count is None else first.start + first.count
    chosen = [a for r in residues[first.start:end] for a in r.atoms.values()]
    return [sum(a.xyz[k] for a in chosen) / len(chosen) for k in range(3)]


def read_session(path: str | Path) -> dict:
    """El JSON, con las estructuras ya leídas (panel["structure"]) y el resaltado con claves
    enteras."""
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if data.get("version") != VERSION or not data.get("panels"):
        raise ValueError(f"Sesión no válida: {path}")
    for panel in data["panels"]:
        panel["structure"] = read_pdb(panel["pdb"])
        panel["highlight"] = {int(k): v for k, v in panel.get("highlight", {}).items()}
    return data
