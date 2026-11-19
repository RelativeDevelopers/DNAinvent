"""Tests de «Cuerpo humano» (dnainvent.anatomy) con un BodyParts3D de juguete: un corazón con
dos piezas, un fémur y una arteria; una tabla ASCT+B pequeña y un Uberon de tres líneas."""

import sqlite3
import zipfile

import numpy as np
import pytest

from dnainvent import anatomy
from dnainvent.anatomy import bodyparts3d, cells, meshes, systems, tissues
from dnainvent.anatomy.build import build_anatomy

CUBE_V = [(0, 0, 0), (10, 0, 0), (10, 10, 0), (0, 10, 0), (0, 0, 10), (10, 0, 10), (10, 10, 10), (0, 10, 10)]
CUBE_F = [(1, 2, 3), (1, 3, 4), (5, 7, 6), (5, 8, 7), (1, 5, 6), (1, 6, 2),
          (2, 6, 7), (2, 7, 3), (3, 7, 8), (3, 8, 4), (4, 8, 5), (4, 5, 1)]


def obj(file_id, fma, name, shift):
    lines = ["#", "# Compatibility version : 4.0", f"# File ID : {file_id}", f"# Concept ID : {fma}",
             f"# English name : {name}",
             f"# Bounds(mm): ({shift},{shift},{shift})-({shift + 10},{shift + 10},{shift + 10})", "#"]
    lines += [f"v {x + shift} {y + shift} {z + shift}" for x, y, z in CUBE_V]
    lines += [f"f {a}//{a} {b}//{b} {c}//{c}" for a, b, c in CUBE_F]
    return "\n".join(lines) + "\n"


def write(path, header, rows):
    path.write_text("\n".join(["\t".join(header)] + ["\t".join(r) for r in rows]) + "\n", encoding="utf-8")


@pytest.fixture
def bp_dir(tmp_path):
    d = tmp_path / "bodyparts3d"
    d.mkdir()
    parts = [("FMA7088", "BP1", "heart"), ("FMA7100", "BP2", "left ventricle"),
             ("FMA9612", "BP3", "left femur"), ("FMA3734", "BP4", "aorta"), ("FMA5022", "BP5", "muscle organ")]
    write(d / "isa_parts_list_e.txt", ["concept id", "representation id", "en"], parts)
    write(d / "partof_parts_list_e.txt", ["concept id", "representation id", "en"], parts[:2])
    write(d / "isa_inclusion_relation_list.txt", ["parent id", "parent name", "child id", "child name"],
          [("FMA5018", "bone organ", "FMA9612", "left femur"), ("FMA50722", "artery", "FMA3734", "aorta")])
    write(d / "partof_inclusion_relation_list.txt", ["parent id", "parent name", "child id", "child name"],
          [("FMA7088", "heart", "FMA7100", "left ventricle"), ("FMA7088", "heart", "FMA9999", "wall of heart")])
    write(d / "isa_element_parts.txt", ["concept id", "name", "element file id"],
          [("FMA7100", "left ventricle", "FJ1"), ("FMA9612", "left femur", "FJ2"), ("FMA3734", "aorta", "FJ3"),
           ("FMA9999", "wall of heart", "FJ4"), ("FMA7088", "heart", "FJ1"), ("FMA7088", "heart", "FJ4")])
    write(d / "partof_element_parts.txt", ["concept id", "name", "element file id"], [("FMA7088", "heart", "FJ1")])
    with zipfile.ZipFile(d / "isa_BP3D_4.0_obj_99.zip", "w") as z:
        z.writestr("isa_BP3D_4.0_obj_99/FJ1.obj", obj("FJ1", "FMA7100", "Left ventricle", 0))
        z.writestr("isa_BP3D_4.0_obj_99/FJ2.obj", obj("FJ2", "FMA9612", "Left femur", 100))
        z.writestr("isa_BP3D_4.0_obj_99/FJ3.obj", obj("FJ3", "FMA3734", "Aorta", 20))
        z.writestr("isa_BP3D_4.0_obj_99/FJ4.obj", obj("FJ4", "", "", 40))        # sin concepto
    (d / "partof_BP3D_4.0_obj_99.zip").write_bytes((d / "isa_BP3D_4.0_obj_99.zip").read_bytes())
    return d


@pytest.fixture
def asctb(tmp_path):
    path = tmp_path / "asct-b-vh-heart.csv"
    path.write_text(
        '"Anatomical Structures, Cell Types and Biomarkers Table for Heart",,,,\n'
        "Version Number:,v1.3,,,\n"
        "AS/1,AS/1/LABEL,AS/1/ID,AS/2,AS/2/LABEL,AS/2/ID,CT/1,CT/1/LABEL,CT/1/ID,BGene/1,BGene/1/LABEL,BGene/1/ID,BGene/2,BGene/2/LABEL,BGene/2/ID\n"
        "heart,heart,UBERON:0000948,,,,cardiac muscle cell,cardiac muscle cell,CL:0000746,TNNT2,troponin T2,HGNC:11949,MYH6,myosin heavy chain 6,HGNC:7576\n"
        "heart,heart,UBERON:0000948,aorta,aorta,FMA:3734,smooth muscle cell,smooth muscle cell,CL:0000192,MYH11,myosin heavy chain 11,HGNC:7569,,,\n",
        encoding="utf-8")
    return cells.read_asctb(path)


@pytest.fixture
def uberon(tmp_path):
    path = tmp_path / "uberon-edit.obo"
    path.write_text("[Term]\nid: UBERON:0000948\nname: heart\nxref: FMA:7088\n\n[Term]\nid: UBERON:0000001\n",
                    encoding="utf-8")
    return cells.read_uberon_fma(path)


def test_bodyparts3d(bp_dir):
    bp = bodyparts3d.BodyParts3D(bp_dir)
    assert set(bp.elements) == {"FJ1", "FJ2", "FJ3", "FJ4"}
    assert bp.elements["FJ1"].fma == "FMA7100" and bp.elements["FJ1"].name == "left ventricle"
    assert bp.elements["FJ4"].fma == "FMA9999"           # el concepto más concreto que la contiene
    assert bp.elements["FJ2"].lo == (100.0, 100.0, 100.0)
    assert bp.partof_parents("FMA7100") == {"FMA7088"} and "FMA5018" in bp.ancestors("FMA9612")
    v, f = bp.mesh("FJ2")
    assert v.shape == (8, 3) and f.shape == (12, 3) and f.min() == 0 and v.min() == 100


def test_missing_files(tmp_path):
    with pytest.raises(FileNotFoundError):
        bodyparts3d.BodyParts3D(tmp_path)


def test_read_obj_splits_polygons():
    v, f = bodyparts3d.read_obj("v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nf 1 2 3 4\n")
    assert f.tolist() == [[0, 1, 2], [0, 2, 3]]


def test_systems(bp_dir):
    bp = bodyparts3d.BodyParts3D(bp_dir)
    assert systems.classify(bp, "FMA9612") == "skeletal"
    assert systems.classify(bp, "FMA3734") == "arteries"
    assert systems.classify(bp, "FMA7100") == "heart"
    assert systems.classify(bp, "FMA9999") == "heart"            # por aquello de lo que forma parte

    class Names:
        def __init__(self, name):
            self.n = name

        def name(self, fma):
            return self.n

        def ancestors(self, fma, tree):
            return set()

        def isa_parents(self, fma):
            return set()

        partof_parents = isa_parents

    for name, expected in [("anterior inferior cerebellar artery", "arteries"), ("left circumflex scapular vein", "veins"),
                           ("flexor retinaculum of right wrist", "connective"), ("left tensor fasciae latae", "muscular"),
                           ("right lens", "senses"), ("left hippocampus", "nervous"), ("gallbladder", "digestive")]:
        assert systems.classify(Names(name), "x") == expected, name


def test_asctb_and_cells(bp_dir, asctb, uberon):
    assert len(asctb) == 2 and asctb[0].structures == [("UBERON:0000948", "heart")]
    assert asctb[0].genes == [("TNNT2", "troponin T2", "HGNC:11949"), ("MYH6", "myosin heavy chain 6", "HGNC:7576")]
    assert asctb[1].structures[-1] == ("FMA3734", "aorta")
    assert uberon == {"FMA7088": {"UBERON:0000948"}}
    bp = bodyparts3d.BodyParts3D(bp_dir)
    index = cells.CellIndex(asctb, uberon)
    found_at, depth, entry = index.find(bp, "FMA7100")          # ventrículo -> corazón (vía Uberon)
    assert (found_at, depth) == ("FMA7088", 1) and "CL:0000746" in entry["cells"]
    assert index.find(bp, "FMA3734")[1] == 0                     # la aorta, directa por FMA
    assert index.find(bp, "FMA9612") is None


def test_tissues():
    assert tissues.tissue_for("muscular", "left biceps") == "skeletal_muscle"
    assert tissues.tissue_for("skeletal", "right sixth costal cartilage") == "cartilage"
    assert tissues.tissue_for("nervous", "left frontal nerve") == "nerve"
    assert tissues.tissue_for("nervous", "left hippocampus") == "brain"
    assert tissues.tissue_for("digestive", "left lower first secondary molar tooth") == "tooth"
    assert tissues.tissue_for("senses", "right lens") is None
    for key, rows in tissues.TISSUES.items():
        assert all(cl.startswith("CL:") and genes for cl, _, genes in rows), key


def test_pack_roundtrip(tmp_path):
    v = np.array(CUBE_V, np.float32) + 5
    f = np.array(CUBE_F) - 1
    q, ff = meshes.weld(meshes.quantize(np.vstack([v, v]), np.zeros(3), 0.5), np.vstack([f, f + 8]))
    assert len(q) == 8 and len(ff) == 24                        # vértices repetidos juntados
    header = meshes.write_pack(tmp_path / "x.dnab", "skeletal", [("FJ9", "FMA1", q, ff)], np.zeros(3), 0.5)
    pack = meshes.read_pack(tmp_path / "x.dnab")
    assert pack.header == header and pack.pieces[0]["id"] == "FJ9"
    verts, faces = pack.mesh(0)
    assert np.allclose(np.sort(verts, axis=0), np.sort(v, axis=0)) and faces.shape == (24, 3)
    normals = meshes.vertex_normals(verts, faces)
    assert np.allclose(np.linalg.norm(normals, axis=1), 1)


def test_split_big_meshes():
    n = meshes.MAX_VERTICES + 10
    q = np.zeros((n, 3), np.uint16)
    faces = np.array([(i, i + 1, i + 2) for i in range(0, n - 2, 3)])
    parts = meshes.split(q, faces)
    assert len(parts) == 2 and sum(len(f) for _, f in parts) == len(faces)
    assert all(len(p) <= meshes.MAX_VERTICES and f.max() < len(p) for p, f in parts)


def test_build_and_open(bp_dir, asctb, uberon, tmp_path, monkeypatch):
    out = tmp_path / "anatomy"
    report = build_anatomy(bodyparts3d.BodyParts3D(bp_dir), asctb, uberon, out, ratio=1.0)
    assert report.pieces == 4 and set(report.systems) == {"heart", "skeletal", "arteries"}
    assert (report.with_cells, report.with_tissue, report.without) == (3, 1, 0)
    monkeypatch.setenv("DNAINVENT_ANATOMY", str(out))
    conn = anatomy.open_anatomy()
    assert conn.execute("SELECT system FROM pieces WHERE id = 'FJ2'").fetchone()[0] == "skeletal"
    row = conn.execute("SELECT s.key, cc.depth FROM concept_cells cc JOIN cell_sets s ON s.id = cc.set_id "
                       "WHERE cc.fma = 'FMA7100'").fetchone()
    assert tuple(row) == ("UBERON:0000948", 1)
    genes = {r[0] for r in conn.execute("SELECT symbol FROM set_cells sc JOIN cell_sets s ON s.id = sc.set_id "
                                        "WHERE s.source = 'asctb' AND s.key = 'UBERON:0000948'")}
    assert genes == {"TNNT2", "MYH6", "MYH11"}      # el corazón entero junta las de sus partes (la aorta)
    kinds = dict(conn.execute("SELECT key, kind FROM cell_sets WHERE source = 'asctb'").fetchall())
    assert kinds == {"UBERON:0000948": "parts", "FMA3734": "leaf"}
    assert conn.execute("SELECT COUNT(*) FROM concept_pieces WHERE fma = 'FMA7088'").fetchone()[0] == 2
    assert conn.execute("SELECT system FROM concepts WHERE fma = 'FMA7088'").fetchone()[0] == "heart"
    assert conn.execute("SELECT name FROM genes WHERE symbol = 'MYH11'").fetchone()[0] == "myosin heavy chain 11"
    assert "CC Attribution 4.0" in conn.execute("SELECT value FROM meta WHERE key = 'bodyparts3d'").fetchone()[0]
    pack = meshes.read_pack(out / "meshes" / "heart.dnab")
    assert {p["id"] for p in pack.pieces} == {"FJ1", "FJ4"}
    with pytest.raises(sqlite3.OperationalError):
        conn.execute("DELETE FROM pieces")                        # solo lectura


def test_open_missing(tmp_path, monkeypatch):
    monkeypatch.setenv("DNAINVENT_ANATOMY", str(tmp_path))
    with pytest.raises(anatomy.AnatomyMissing):
        anatomy.open_anatomy()
