"""Lectura de BodyParts3D 4.0: el cuerpo de referencia (hombre adulto) pieza a pieza.

    bp = BodyParts3D("C:/.../datos/bodyparts3d")
    bp.elements["FJ1252"]          # Element(file_id='FJ1252', fma='FMA59763', name='Gingiva of upper jaw', ...)
    bp.name("FMA7088")             # 'heart'
    bp.partof_parents("FMA59763")  # conceptos de los que forma parte
    vertices, faces = bp.mesh("FJ1252")

La carpeta tiene los archivos tal como se descargan de
https://dbarchive.biosciencedbc.jp/data/bodyparts3d/LATEST/ (no hace falta descomprimir):

    isa_parts_list_e.txt, partof_parts_list_e.txt          concepto (FMA) -> nombre en inglés
    isa_inclusion_relation_list.txt, partof_...             padre -> hijo (dos árboles: «es un»
                                                            y «forma parte de»)
    isa_element_parts.txt, partof_element_parts.txt         concepto -> piezas 3D que lo forman
    isa_BP3D_4.0_obj_99.zip                                 una malla OBJ por pieza (2.234), en mm

Cada pieza («element», archivo FJxxxx.obj) es un trozo que se puede seleccionar. Los conceptos
pueden juntar varias piezas (el corazón, «hueso», «arteria»...). Los identificadores son los de
la Foundational Model of Anatomy (FMA), que es lo que permite enlazar con otras bases de datos.

Licencia: CC BY 4.0 desde el 27-02-2025 (la cabecera de los OBJ todavía dice la antigua, CC BY-SA
2.1 Japón; manda la página de la licencia). Hay que citar siempre ATTRIBUTION.
"""

import csv
import re
import zipfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

ATTRIBUTION = ("BodyParts3D, © The Database Center for Life Science licensed under CC Attribution "
               "4.0 International")
LICENSE_URL = "https://dbarchive.biosciencedbc.jp/en/bodyparts3d/lic.html"
MESH_ZIP = "isa_BP3D_4.0_obj_99.zip"
TABLES = ("isa_parts_list_e.txt", "partof_parts_list_e.txt", "isa_inclusion_relation_list.txt",
          "partof_inclusion_relation_list.txt", "isa_element_parts.txt", "partof_element_parts.txt")

_HEADER = re.compile(r"^# (File ID|Concept ID|English name|Bounds\(mm\)) ?: ?(.*)$", re.M)


@dataclass
class Element:
    file_id: str
    fma: str                 # el concepto de la pieza (FMA...), deducido si el OBJ no lo trae
    name: str                # nombre en inglés
    lo: tuple                # caja que la contiene, en mm
    hi: tuple


def read_table(path: Path) -> list[list[str]]:
    """Las tablas son texto separado por tabuladores con una fila de cabecera."""
    with open(path, encoding="utf-8", newline="") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    return [r for r in rows[1:] if r]


def parse_header(text: str) -> dict:
    """La cabecera de comentarios de un OBJ de BodyParts3D."""
    out = {}
    for key, value in _HEADER.findall(text[:3000]):
        out[key] = value.strip()
    if "Bounds(mm)" in out:
        nums = [float(x) for x in re.findall(r"-?\d+(?:\.\d+)?", out["Bounds(mm)"])]
        out["lo"], out["hi"] = tuple(nums[:3]), tuple(nums[3:6])
    return out


def read_obj(text: str):
    """Vértices (N×3, float32, mm) y triángulos (M×3, int32, desde 0) de un OBJ.
    Se ignoran las normales: el visor las calcula."""
    import numpy as np
    verts, faces = [], []
    for line in text.splitlines():
        if line.startswith("v "):
            verts.append(line.split()[1:4])
        elif line.startswith("f "):
            idx = [int(p.split("/")[0]) - 1 for p in line.split()[1:]]
            for k in range(1, len(idx) - 1):            # los polígonos se parten en triángulos
                faces.append((idx[0], idx[k], idx[k + 1]))
    return (np.asarray(verts, dtype=np.float32).reshape(-1, 3),
            np.asarray(faces, dtype=np.int32).reshape(-1, 3))


class BodyParts3D:
    def __init__(self, folder):
        self.folder = Path(folder)
        missing = [n for n in TABLES + (MESH_ZIP,) if not (self.folder / n).is_file()]
        if missing:
            raise FileNotFoundError(f"Faltan archivos de BodyParts3D en {self.folder}: {', '.join(missing)}")
        self._names: dict[str, str] = {}
        for tree in ("isa", "partof"):
            for row in read_table(self.folder / f"{tree}_parts_list_e.txt"):
                self._names.setdefault(row[0], row[2])
        self._parents = {"isa": defaultdict(set), "partof": defaultdict(set)}
        for tree in ("isa", "partof"):
            for parent, pname, child, cname in read_table(self.folder / f"{tree}_inclusion_relation_list.txt"):
                self._parents[tree][child].add(parent)
                self._names.setdefault(parent, pname)
                self._names.setdefault(child, cname)
        self.members: dict[str, set[str]] = defaultdict(set)      # concepto -> piezas
        for tree in ("isa", "partof"):
            for concept, name, file_id in read_table(self.folder / f"{tree}_element_parts.txt"):
                self.members[concept].add(file_id)
                self._names.setdefault(concept, name)
        self._zip = zipfile.ZipFile(self.folder / MESH_ZIP)
        self._paths = {Path(i.filename).stem: i.filename for i in self._zip.infolist()
                       if i.filename.lower().endswith(".obj")}
        self.elements = self._read_elements()

    # ---- nombres y árboles
    def name(self, fma: str) -> str:
        return self._names.get(fma, "")

    def concepts(self) -> dict[str, str]:
        return dict(self._names)

    def isa_parents(self, fma: str) -> set[str]:
        return set(self._parents["isa"].get(fma, ()))

    def partof_parents(self, fma: str) -> set[str]:
        return set(self._parents["partof"].get(fma, ()))

    def ancestors(self, fma: str, tree: str = "isa") -> set[str]:
        out, stack = set(), [fma]
        while stack:
            for parent in self._parents[tree].get(stack.pop(), ()):
                if parent not in out:
                    out.add(parent)
                    stack.append(parent)
        return out

    def relations(self, tree: str):
        """(hijo, padre) de un árbol, para guardarlo."""
        for child, parents in self._parents[tree].items():
            for parent in parents:
                yield child, parent

    # ---- piezas
    def _read_elements(self) -> dict[str, Element]:
        elements = {}
        for file_id, path in sorted(self._paths.items()):
            with self._zip.open(path) as f:
                head = parse_header(f.read(3000).decode("utf-8", "replace"))
            fma = head.get("Concept ID", "")
            if not fma:
                # 17 piezas vienen sin concepto: el más concreto que la contiene
                owners = [c for c, files in self.members.items() if file_id in files]
                fma = min(owners, key=lambda c: (len(self.members[c]), c)) if owners else ""
            name = self._names.get(fma) or head.get("English name", "")
            elements[file_id] = Element(file_id, fma, name, head.get("lo", (0, 0, 0)), head.get("hi", (0, 0, 0)))
        return elements

    def mesh(self, file_id: str):
        with self._zip.open(self._paths[file_id]) as f:
            return read_obj(f.read().decode("utf-8", "replace"))
