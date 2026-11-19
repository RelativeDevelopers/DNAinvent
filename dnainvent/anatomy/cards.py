"""Las fichas de «Cuerpo humano»: buscar una parte del cuerpo y ver qué células tiene y qué genes
las marcan. Lo usan igual la consola y la web.

    conn = anatomy.open_anatomy()
    cards.search(conn, "corazón", "es")                 # [{'fma': 'FMA7088', ...}, ...]
    card = cards.concept_card(conn, "FMA7088", "es", series={"HBB", "DMD", ...})
    cards.piece_concept(conn, "FJ3259")                 # 'FMA24475' (la pieza 3D -> su concepto)
    cards.concept_pieces(conn, "FMA7088")               # las 83 piezas del corazón

Un «concepto» es una parte del cuerpo con nombre (el corazón, el fémur izquierdo, el hueso...)
y puede juntar muchas piezas 3D. Las células salen, por orden:
    own     las tablas ASCT+B tienen esa parte en concreto;
    organ   las tienen para algo de lo que forma parte (el ventrículo izquierdo -> el corazón);
    tissue  el tejido general de su sistema (tissues.py: «células típicas del hueso»);
    None    todavía no hay datos (los sentidos, el aparato digestivo...).
"""

from collections import Counter

from . import names
from .systems import SYSTEM_KEYS
from .tissues import CL_URL, SOURCES as TISSUE_SOURCES

BODYPARTS3D_URL = "https://dbarchive.biosciencedbc.jp/en/bodyparts3d/"
HRA_URL = "https://humanatlas.io/"
UBERON_URL = "https://www.ebi.ac.uk/ols4/ontologies/uberon/classes?obo_id={id}"
NCBI_GENE_URL = "https://www.ncbi.nlm.nih.gov/gene/?term={symbol}%5Bsym%5D+AND+human%5Borgn%5D"
PARTS_SHOWN = 40
BROAD = 400            # un concepto con más piezas que esto es una clase («hueso», «arteria»)


# ---------------------------------------------------------------- sistemas y búsqueda
def systems(conn) -> list[dict]:
    return [dict(r) for r in conn.execute("SELECT * FROM systems ORDER BY ord")]


def meta(conn) -> dict:
    return {r["key"]: r["value"] for r in conn.execute("SELECT key, value FROM meta")}


_counts: dict[int, dict] = {}


def _piece_counts(conn) -> dict[str, int]:
    """Cuántas piezas tiene cada concepto (una vez por conexión)."""
    key = id(conn)
    if key not in _counts:
        _counts[key] = {r[0]: r[1] for r in conn.execute("SELECT fma, COUNT(*) FROM concept_pieces GROUP BY fma")}
    return _counts[key]


def concept_name(conn, fma: str) -> str:
    row = conn.execute("SELECT name FROM concepts WHERE fma = ?", (fma,)).fetchone()
    return row[0] if row else fma


def label(conn, fma: str, language: str) -> dict:
    """{'fma', 'name' (inglés), 'label' (escrito a mano en el idioma o None)}."""
    return {"fma": fma, "name": names.english(concept_name(conn, fma)),
            "label": names.concept(fma, language) if language != "en" else None}


def _main_system(conn, fma: str) -> str | None:
    """El sistema de la parte del cuerpo (systems.classify al construir); si no encaja en
    ninguno, el de la mayoría de sus piezas."""
    own = conn.execute("SELECT system FROM concepts WHERE fma = ?", (fma,)).fetchone()
    if own and own[0] and own[0] != "other":
        return own[0]
    row = conn.execute("SELECT p.system, COUNT(*) n FROM concept_pieces cp JOIN pieces p ON p.id = cp.piece "
                       "WHERE cp.fma = ? GROUP BY p.system ORDER BY n DESC LIMIT 1", (fma,)).fetchone()
    return row[0] if row else None


def search(conn, text: str, language: str, limit: int = 25) -> list[dict]:
    """Partes del cuerpo (con alguna pieza 3D) cuyo nombre en el idioma o en inglés contiene
    todas las palabras buscadas. Primero las que empiezan igual y las que tienen nombre propio
    en el idioma; las clases muy generales («hueso»), al final."""
    words = names.normalize(text).split()
    if not words or len("".join(words)) < 2:
        return []
    counts = _piece_counts(conn)
    own = names.concepts(language) if language != "en" else {}
    found = []
    for fma, name in conn.execute("SELECT fma, name FROM concepts"):
        if fma not in counts:
            continue
        candidates = [names.normalize(own[fma])] if fma in own else []
        candidates.append(names.normalize(name))
        best = None
        for c in candidates:
            if not all(w in c for w in words):
                continue
            joined = " ".join(words)
            score = 0 if c == joined else 1 if c.startswith(joined) else 2 if f" {words[0]}" in f" {c}" else 3
            best = score if best is None else min(best, score)
        if best is not None:
            found.append((best, fma not in own, counts[fma] > BROAD, len(name), name, fma))
    found.sort()
    out = []
    for *_, fma in found[:limit]:
        item = label(conn, fma, language)
        item.update(pieces=counts[fma], system=_main_system(conn, fma))
        out.append(item)
    return out


def concept_pieces(conn, fma: str) -> list[str]:
    return [r[0] for r in conn.execute("SELECT piece FROM concept_pieces WHERE fma = ? ORDER BY piece", (fma,))]


def piece_concept(conn, piece_id: str) -> str | None:
    row = conn.execute("SELECT fma FROM pieces WHERE id = ?", (piece_id,)).fetchone()
    return row[0] if row else None


def system_concepts(conn, system: str, language: str, text: str = "") -> list[dict]:
    """Las piezas de un sistema, por su concepto (una fila por concepto), en orden alfabético
    del nombre que se ve; text filtra como search."""
    words = names.normalize(text).split()
    rows = conn.execute("SELECT DISTINCT fma FROM pieces WHERE system = ?", (system,)).fetchall()
    out = []
    for (fma,) in rows:
        item = label(conn, fma, language)
        shown = names.normalize(item["label"] or item["name"])
        if words and not all(w in shown or w in names.normalize(item["name"]) for w in words):
            continue
        out.append(item)
    out.sort(key=lambda i: names.normalize(i["label"] or i["name"]))
    return out


# ---------------------------------------------------------------- la ficha
def _relatives(conn, fma: str, column: str, other: str, language: str) -> list[dict]:
    """Padres (de qué forma parte) o hijos (qué partes tiene) en el árbol «forma parte de», solo
    los que tienen piezas 3D."""
    counts = _piece_counts(conn)
    rows = conn.execute(f"SELECT DISTINCT {other} FROM relations WHERE tree = 'partof' AND {column} = ?", (fma,))
    out = [label(conn, r[0], language) for r in rows if r[0] in counts and r[0] != fma]
    return sorted(out, key=lambda i: names.normalize(i["label"] or i["name"]))


def kinds(conn, fma: str, language: str) -> list[dict]:
    """Qué es («es un»): las clases inmediatas del árbol «es un»."""
    rows = conn.execute("SELECT parent FROM relations WHERE tree = 'isa' AND child = ? ORDER BY parent", (fma,))
    return [label(conn, r[0], language) for r in rows]


def _cell_set(conn, set_id: int, language: str, series) -> list[dict]:
    cells: dict[str, dict] = {}
    rows = conn.execute("SELECT sc.cl, ct.name, sc.symbol, g.name FROM set_cells sc "
                        "LEFT JOIN cell_types ct ON ct.cl = sc.cl LEFT JOIN genes g ON g.symbol = sc.symbol "
                        "WHERE sc.set_id = ? ORDER BY sc.rowid", (set_id,))
    for cl, cell_name, symbol, gene_name in rows:
        cell = cells.setdefault(cl, {"cl": cl, "name": cell_name or cl,
                                     "label": names.cell(cl, language) if language != "en" else None,
                                     "url": CL_URL.format(cl=cl) if cl.startswith("CL:") else None,
                                     "genes": []})
        if symbol and all(g["symbol"] != symbol for g in cell["genes"]):
            cell["genes"].append({"symbol": symbol, "name": gene_name, "series": symbol in series})
    # primero las que tienen más genes marcadores: suelen ser las propias del órgano (el
    # hepatocito en el hígado), y las de todas partes (fibroblastos, vasos) al final
    return sorted(cells.values(), key=lambda c: (-len(c["genes"]), c["name"].lower()))


def cells_for(conn, fma: str, language: str, series=frozenset()) -> dict:
    """Las células de una parte del cuerpo y de dónde salen (ver el principio del módulo)."""
    row = conn.execute("SELECT cc.found_at, cc.depth, s.id, s.key, s.label, s.organ, s.kind FROM concept_cells cc "
                       "JOIN cell_sets s ON s.id = cc.set_id WHERE cc.fma = ?", (fma,)).fetchone()
    if row is not None:
        found_at, depth, set_id, key, set_label, organ, kind = row
        return {"source": "own" if found_at == fma else "organ", "in_parts": kind == "parts",
                "found_at": label(conn, found_at, language), "depth": depth, "key": key,
                "uberon_url": UBERON_URL.format(id=key) if key.startswith("UBERON:") else None,
                "table": organ, "cells": _cell_set(conn, set_id, language, series)}
    tissue = conn.execute("SELECT p.tissue, COUNT(*) n FROM concept_pieces cp JOIN pieces p ON p.id = cp.piece "
                          "WHERE cp.fma = ? AND p.tissue IS NOT NULL GROUP BY p.tissue ORDER BY n DESC LIMIT 1",
                          (fma,)).fetchone()
    if tissue:
        set_row = conn.execute("SELECT id FROM cell_sets WHERE source = 'tissue' AND key = ?", (tissue[0],)).fetchone()
        if set_row:
            return {"source": "tissue", "tissue": tissue[0], "cells": _cell_set(conn, set_row[0], language, series)}
    return {"source": None, "cells": []}


def concept_card(conn, fma: str, language: str, series=frozenset()) -> dict | None:
    """Todo lo de una parte del cuerpo. series: los símbolos de los genes de serie (para marcar
    los que tienen ficha completa en «Genes reales»)."""
    if not conn.execute("SELECT 1 FROM concepts WHERE fma = ?", (fma,)).fetchone():
        return None
    pieces = conn.execute("SELECT p.id, p.system, p.lo_x, p.lo_y, p.lo_z, p.hi_x, p.hi_y, p.hi_z FROM concept_pieces cp "
                          "JOIN pieces p ON p.id = cp.piece WHERE cp.fma = ?", (fma,)).fetchall()
    card = label(conn, fma, language)
    systems_count = Counter(p["system"] for p in pieces)
    card["pieces"] = len(pieces)
    card["systems"] = [k for k in SYSTEM_KEYS if k in systems_count]
    card["system"] = _main_system(conn, fma) if pieces else None
    card["size_cm"] = size_cm(pieces)
    card["kinds"] = kinds(conn, fma, language)
    card["part_of"] = _relatives(conn, fma, "child", "parent", language)
    parts = _relatives(conn, fma, "parent", "child", language)
    card["parts"], card["more_parts"] = parts[:PARTS_SHOWN], max(0, len(parts) - PARTS_SHOWN)
    card["cells"] = cells_for(conn, fma, language, series)
    genes = {}
    for cell in card["cells"]["cells"]:
        for g in cell["genes"]:
            genes.setdefault(g["symbol"], g)
    card["genes"] = list(genes.values())
    card["series_genes"] = [g["symbol"] for g in card["genes"] if g["series"]]
    card["links"] = links(card)
    return card


def size_cm(pieces) -> float | None:
    """Lo que mide (el lado más largo de su caja), en cm. Si son una pieza de cada lado del
    cuerpo (los dos fémures), lo que mide una de ellas, no la distancia entre las dos."""
    if not pieces:
        return None
    left = [p for p in pieces if p["hi_x"] < 0]
    right = [p for p in pieces if p["lo_x"] > 0]
    if left and right and len(left) + len(right) == len(pieces):
        pieces = max(left, right, key=len)

    def extent(group):
        lo = [min(p[k] for p in group) for k in ("lo_x", "lo_y", "lo_z")]
        hi = [max(p[k] for p in group) for k in ("hi_x", "hi_y", "hi_z")]
        return max(h - l for h, l in zip(hi, lo))
    return round(extent(pieces) / 10, 1)


def links(card: dict) -> list[tuple[str, str]]:
    out = [("BodyParts3D", BODYPARTS3D_URL)]
    cells = card["cells"]
    if cells.get("uberon_url"):
        out.append((f"Uberon {cells['key']}", cells["uberon_url"]))
    if cells["source"] in ("own", "organ"):
        out.append(("Human Reference Atlas", HRA_URL))
    elif cells["source"] == "tissue":
        out += [("Cell Ontology", TISSUE_SOURCES[0]), ("Human Protein Atlas", TISSUE_SOURCES[1])]
    return out


def gene_url(symbol: str) -> str:
    return NCBI_GENE_URL.format(symbol=symbol)
