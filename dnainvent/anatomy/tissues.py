"""El tejido general de cada sistema, para las piezas que no tienen células propias en las
tablas ASCT+B (muchos músculos, huesos, vasos y nervios concretos).

Es una tabla escrita a mano, como curated.py en la base de datos de genes: pocos tipos de
célula por tejido, los principales, y genes marcadores clásicos (los que se usan para
reconocerlos en los libros y en los atlas de células). La ficha lo dice claro: «células típicas
del músculo esquelético», no «células de este músculo en concreto».

Identificadores: Cell Ontology (CL) para las células y símbolos de HGNC para los genes.
Fuentes: la Cell Ontology (definición de cada tipo) y el atlas de tipos de célula del Human
Protein Atlas (genes marcadores).
"""

import re

CL_URL = "https://www.ebi.ac.uk/ols4/ontologies/cl/classes?obo_id={cl}"
SOURCES = ("https://obofoundry.org/ontology/cl.html",
           "https://www.proteinatlas.org/humanproteome/single+cell")

# tejido -> [(CL, nombre en inglés, [genes])]
TISSUES = {
    "skeletal_muscle": [
        ("CL:0008002", "skeletal muscle fiber", ["ACTA1", "DES", "TTN", "DMD", "MYH7", "MYH2"]),
        ("CL:0000594", "skeletal muscle satellite cell", ["PAX7", "MYF5"]),
        ("CL:0000057", "fibroblast", ["COL1A1", "PDGFRA", "DCN"]),
    ],
    "bone": [
        ("CL:0000062", "osteoblast", ["RUNX2", "SP7", "COL1A1", "BGLAP"]),
        ("CL:0000137", "osteocyte", ["SOST", "DMP1"]),
        ("CL:0000092", "osteoclast", ["ACP5", "CTSK"]),
    ],
    "cartilage": [
        ("CL:0000138", "chondrocyte", ["SOX9", "COL2A1", "ACAN"]),
    ],
    "tooth": [
        ("CL:0000060", "odontoblast", ["DSPP", "DMP1"]),
    ],
    "dense_connective": [
        ("CL:0000388", "tendon cell", ["SCX", "TNMD"]),
        ("CL:0000057", "fibroblast", ["COL1A1", "PDGFRA", "DCN"]),
    ],
    "artery": [
        ("CL:1000413", "endothelial cell of artery", ["PECAM1", "VWF", "CDH5"]),
        ("CL:0000192", "smooth muscle cell", ["MYH11", "ACTA2", "TAGLN"]),
    ],
    "vein": [
        ("CL:0002543", "vein endothelial cell", ["PECAM1", "VWF", "ACKR1"]),
        ("CL:0000192", "smooth muscle cell", ["MYH11", "ACTA2", "TAGLN"]),
    ],
    "nerve": [
        ("CL:0002573", "Schwann cell", ["MPZ", "PMP22", "S100B"]),
    ],
    "brain": [
        ("CL:0000540", "neuron", ["RBFOX3", "SNAP25"]),
        ("CL:0000127", "astrocyte", ["GFAP", "AQP4"]),
        ("CL:0000128", "oligodendrocyte", ["MBP", "OLIG2", "MOG"]),
        ("CL:0000129", "microglial cell", ["P2RY12", "AIF1"]),
    ],
    "skin": [
        ("CL:0000312", "keratinocyte", ["KRT5", "KRT14", "KRT10"]),
        ("CL:0000148", "melanocyte", ["PMEL", "MLANA", "TYR"]),
    ],
}

_NERVE = re.compile(r"nerve|ganglion|ramus|plexus of nerves")
_CARTILAGE = re.compile(r"cartilage|meniscus|intervertebral|articular disk")
_TOOTH = re.compile(r"tooth")


def tissue_for(system: str, name: str) -> str | None:
    """El tejido general de una pieza según su sistema y su nombre (en inglés), o None si no
    hay uno que tenga sentido (los sentidos, el aparato digestivo...: ahí se dice que aún no hay
    datos)."""
    name = name.lower()
    if _TOOTH.search(name):
        return "tooth"
    if system == "muscular":
        return "skeletal_muscle"
    if system == "skeletal":
        return "cartilage" if _CARTILAGE.search(name) else "bone"
    if system == "connective":
        return "dense_connective"
    if system == "arteries":
        return "artery"
    if system == "veins":
        return "vein"
    if system == "nervous":
        return "nerve" if _NERVE.search(name) else "brain"
    if system == "skin":
        return "skin"
    return None
