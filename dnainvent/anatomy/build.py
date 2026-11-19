"""Construye «Cuerpo humano»: la base de datos y las mallas empaquetadas.

    report = build_anatomy(BodyParts3D(carpeta), asctb_rows, uberon_fma, out_dir)

Resultado (out_dir, por defecto ~/.dnainvent/anatomy):
    anatomy.db          piezas, conceptos y sus árboles, sistemas, células y genes
    meshes/<sistema>.dnab   las mallas de cada sistema (formato en meshes.py)

Se hace en el ordenador de quien publica (tools/build_anatomy.py), como el catálogo de genes:
la consola y la web solo lo leen.
"""

import json
import sqlite3
import time
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np

from .. import __version__
from . import FORMAT_VERSION
from .bodyparts3d import ATTRIBUTION, LICENSE_URL
from .cells import ASCTB_CITATION, CellIndex
from .meshes import quantize, simplify, weld, write_pack
from .systems import SYSTEMS, classify
from .tissues import SOURCES as TISSUE_SOURCES, TISSUES, tissue_for

DEFAULT_RATIO = 0.25          # se queda con ~1/4 de los triángulos (ver Decisiones en la bitácora)
STEP_MM = 0.03                # resolución de las posiciones guardadas

SCHEMA = """
CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT);
CREATE TABLE systems (key TEXT PRIMARY KEY, ord INTEGER, color TEXT, pieces INTEGER,
                      triangles INTEGER, file TEXT, bytes INTEGER);
CREATE TABLE concepts (fma TEXT PRIMARY KEY, name TEXT NOT NULL, system TEXT);
CREATE TABLE relations (tree TEXT NOT NULL, child TEXT NOT NULL, parent TEXT NOT NULL,
                        PRIMARY KEY (tree, child, parent));
CREATE INDEX relations_parent ON relations (parent);
CREATE TABLE concept_pieces (fma TEXT NOT NULL, piece TEXT NOT NULL, PRIMARY KEY (fma, piece));
CREATE INDEX concept_pieces_piece ON concept_pieces (piece);
CREATE TABLE pieces (id TEXT PRIMARY KEY, fma TEXT NOT NULL, system TEXT NOT NULL, tissue TEXT,
                     triangles INTEGER, lo_x REAL, lo_y REAL, lo_z REAL, hi_x REAL, hi_y REAL, hi_z REAL);
CREATE INDEX pieces_fma ON pieces (fma);
CREATE TABLE cell_sets (id INTEGER PRIMARY KEY, source TEXT NOT NULL, key TEXT NOT NULL,
                        label TEXT, organ TEXT, kind TEXT, UNIQUE (source, key, kind));
CREATE TABLE cell_types (cl TEXT PRIMARY KEY, name TEXT NOT NULL);
CREATE TABLE genes (symbol TEXT PRIMARY KEY, name TEXT, hgnc TEXT);
CREATE TABLE set_cells (set_id INTEGER NOT NULL, cl TEXT NOT NULL, symbol TEXT);
CREATE INDEX set_cells_set ON set_cells (set_id);
CREATE TABLE concept_cells (fma TEXT PRIMARY KEY, set_id INTEGER NOT NULL, found_at TEXT,
                            depth INTEGER);
"""

Progress = Callable[[str, dict], None]


@dataclass
class BuildReport:
    pieces: int = 0
    triangles: int = 0
    bytes: int = 0
    systems: dict = field(default_factory=dict)
    with_cells: int = 0           # piezas con células propias de ASCT+B (o de lo que forman parte)
    with_tissue: int = 0          # piezas que se quedan con el tejido general
    without: int = 0
    seconds: float = 0.0


def build_anatomy(bp, asctb_rows, uberon_fma: dict, out_dir, ratio: float = DEFAULT_RATIO,
                  progress: Progress | None = None) -> BuildReport:
    started = time.time()
    progress = progress or (lambda key, data: None)
    out_dir = Path(out_dir)
    (out_dir / "meshes").mkdir(parents=True, exist_ok=True)
    report = BuildReport()

    # ---- piezas y sistemas
    pieces = []
    for element in sorted(bp.elements.values(), key=lambda e: e.file_id):
        system = classify(bp, element.fma)
        pieces.append((element, system, tissue_for(system, element.name)))
    lo = np.min([e.lo for e, _, _ in pieces], axis=0) - 1.0
    origin = np.floor(lo)

    # ---- mallas
    by_system = defaultdict(list)
    triangles = {}
    for n, (element, system, _) in enumerate(pieces, 1):
        vertices, faces = bp.mesh(element.file_id)
        vertices, faces = simplify(vertices, faces, ratio)
        q, f = weld(quantize(vertices, origin, STEP_MM), faces)
        by_system[system].append((element.file_id, element.fma, q, f))
        triangles[element.file_id] = len(f)
        if n % 50 == 0 or n == len(pieces):
            progress("anatomy.build.meshes", {"done": n, "total": len(pieces)})
    systems_rows = []
    for order, (key, color) in enumerate(SYSTEMS):
        items = by_system.get(key)
        if not items:
            continue
        path = out_dir / "meshes" / f"{key}.dnab"
        write_pack(path, key, items, origin, STEP_MM)
        size = path.stat().st_size
        tris = sum(len(f) for _, _, _, f in items)
        systems_rows.append((key, order, color, len(items), tris, f"meshes/{key}.dnab", size))
        report.systems[key] = {"pieces": len(items), "triangles": tris, "bytes": size}
        report.bytes += size
        report.triangles += tris
    report.pieces = len(pieces)

    # ---- células
    index = CellIndex(asctb_rows, uberon_fma)
    db_path = out_dir / "anatomy.db"
    tmp = db_path.with_suffix(".part")
    tmp.unlink(missing_ok=True)
    conn = sqlite3.connect(tmp)
    conn.executescript(SCHEMA)
    with conn:
        meta = {"format": str(FORMAT_VERSION), "dnainvent": __version__, "built": time.strftime("%Y-%m-%d"),
                "ratio": str(ratio), "step_mm": str(STEP_MM), "origin": json.dumps(origin.tolist()),
                "bodyparts3d": ATTRIBUTION, "bodyparts3d_license": LICENSE_URL, "asctb": ASCTB_CITATION,
                "tissue_sources": json.dumps(TISSUE_SOURCES)}
        conn.executemany("INSERT INTO meta VALUES (?, ?)", meta.items())
        conn.executemany("INSERT INTO systems VALUES (?, ?, ?, ?, ?, ?, ?)", systems_rows)
        # las piezas de cada concepto (el corazón: 83), según las dos tablas «element_parts»,
        # más la propia de cada pieza
        members = {(concept, file_id) for concept, files in bp.members.items() for file_id in files
                   if file_id in bp.elements}
        members |= {(e.fma, e.file_id) for e in bp.elements.values()}
        with_pieces = {concept for concept, _ in members}
        conn.executemany("INSERT INTO concepts VALUES (?, ?, ?)",
                         [(fma, name, classify(bp, fma) if fma in with_pieces else None)
                          for fma, name in sorted(bp.concepts().items())])
        for tree in ("isa", "partof"):
            conn.executemany("INSERT OR IGNORE INTO relations VALUES (?, ?, ?)",
                             ((tree, c, p) for c, p in bp.relations(tree)))
        conn.executemany("INSERT INTO concept_pieces VALUES (?, ?)", sorted(members))
        conn.executemany("INSERT INTO pieces VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                         [(e.file_id, e.fma, system, tissue, triangles[e.file_id], *e.lo, *e.hi)
                          for e, system, tissue in pieces])

        cell_names, genes, set_ids = {}, {}, {}

        def cell_set(source, key, label, organ, cells, kind=None):
            if (source, key, kind) in set_ids:
                return set_ids[(source, key, kind)]
            cur = conn.execute("INSERT INTO cell_sets (source, key, label, organ, kind) VALUES (?, ?, ?, ?, ?)",
                               (source, key, label, organ, kind))
            set_id = set_ids[(source, key, kind)] = cur.lastrowid
            for cl, name, markers in cells:
                cell_names.setdefault(cl, name)
                rows = [(set_id, cl, symbol) for symbol in markers] or [(set_id, cl, None)]
                conn.executemany("INSERT INTO set_cells VALUES (?, ?, ?)", rows)
            return set_id

        for key, cells in TISSUES.items():
            for _, _, markers in cells:
                for symbol in markers:
                    genes.setdefault(symbol, (None, None))
            cell_set("tissue", key, key, None, cells)

        # células de cada concepto con piezas (las piezas y lo que las junta: el corazón entero
        # también, no solo sus trozos)
        for fma in sorted(with_pieces):
            hit = index.find(bp, fma)
            if hit:
                found_at, depth, entry = hit
                cells = []
                for cl, cell in sorted(entry["cells"].items()):
                    for symbol, (name, hgnc) in cell["genes"].items():
                        if genes.get(symbol, (None, None))[0] is None:
                            genes[symbol] = (name, hgnc)
                    cells.append((cl, cell["label"], sorted(cell["genes"])))
                set_id = cell_set("asctb", entry["key"], entry["label"], entry["organ"], cells, entry["kind"])
                conn.execute("INSERT INTO concept_cells VALUES (?, ?, ?, ?)", (fma, set_id, found_at, depth))
        conn.executemany("INSERT INTO cell_types VALUES (?, ?)", sorted(cell_names.items()))
        conn.executemany("INSERT INTO genes VALUES (?, ?, ?)",
                         sorted((s, n, h) for s, (n, h) in genes.items()))
    with_cells = {r[0] for r in conn.execute("SELECT fma FROM concept_cells")}
    conn.close()
    tmp.replace(db_path)

    for element, _, tissue in pieces:
        if element.fma in with_cells:
            report.with_cells += 1
        elif tissue:
            report.with_tissue += 1
        else:
            report.without += 1
    report.bytes += db_path.stat().st_size
    report.seconds = round(time.time() - started, 1)
    progress("anatomy.build.done", {"pieces": report.pieces})
    return report
