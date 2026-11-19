"""Construye «Cuerpo humano»: el cuerpo 3D pieza a pieza, con sus células y sus genes.

    python tools/build_anatomy.py                           # -> ~/.dnainvent/anatomy/
    python tools/build_anatomy.py --bodyparts3d C:/ruta/bodyparts3d
    python tools/build_anatomy.py --ratio 0.5               # más detalle (archivos más grandes)
    python tools/build_anatomy.py --out C:/ruta/anatomy

Necesita:
  * Los 8 archivos de BodyParts3D 4.0 (https://dbarchive.biosciencedbc.jp/data/bodyparts3d/LATEST/):
    isa_/partof_parts_list_e.txt, isa_/partof_inclusion_relation_list.txt,
    isa_/partof_element_parts.txt e isa_/partof_BP3D_4.0_obj_99.zip, sin descomprimir. Por
    defecto se buscan en ../datos/bodyparts3d (al lado de la carpeta del programa).
  * Conexión, la primera vez: las tablas ASCT+B del Human Reference Atlas y Uberon se descargan
    de GitHub y se guardan en ~/.dnainvent/cache/anatomy/.
  * pip install numpy fast-simplification   (para reducir los triángulos)
"""

import argparse
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dnainvent.anatomy import anatomy_dir  # noqa: E402
from dnainvent.anatomy.bodyparts3d import ATTRIBUTION, BodyParts3D  # noqa: E402
from dnainvent.anatomy.build import DEFAULT_RATIO, build_anatomy  # noqa: E402
from dnainvent.anatomy.cells import (ASCTB_CITATION, ASCTB_ORGANS, ASCTB_URL, UBERON_URL,  # noqa: E402
                                     read_asctb, read_uberon_fma)

CACHE = Path.home() / ".dnainvent" / "cache" / "anatomy"
AGENT = "DNAinvent (Relative Developers; https://github.com/RelativeDevelopers/DNAinvent)"


def download(url: str, target: Path) -> Path:
    """Descarga una vez y guarda; las siguientes veces usa la copia."""
    if target.is_file() and target.stat().st_size > 0:
        return target
    target.parent.mkdir(parents=True, exist_ok=True)
    print(f"  descargando {url.rsplit('/', 1)[-1]}...", flush=True)
    request = urllib.request.Request(url, headers={"User-Agent": AGENT})
    with urllib.request.urlopen(request, timeout=120) as response:
        data = response.read()
    partial = target.with_suffix(target.suffix + ".part")
    partial.write_bytes(data)
    partial.replace(target)
    time.sleep(0.3)                    # sin prisas: son pocos archivos
    return target


def main() -> int:
    parser = argparse.ArgumentParser(description="Construye «Cuerpo humano» (anatomía 3D, células y genes).")
    parser.add_argument("--bodyparts3d", default=str(ROOT.parent / "datos" / "bodyparts3d"),
                        help="carpeta con los archivos de BodyParts3D 4.0")
    parser.add_argument("--out", default=str(anatomy_dir()), help="carpeta de salida (por defecto ~/.dnainvent/anatomy)")
    parser.add_argument("--ratio", type=float, default=DEFAULT_RATIO,
                        help=f"parte de los triángulos que se conserva (por defecto {DEFAULT_RATIO})")
    args = parser.parse_args()

    try:
        import fast_simplification  # noqa: F401
    except ImportError:
        print("Falta fast-simplification: pip install fast-simplification")
        return 1

    print(f"Cuerpo humano de DNAinvent -> {args.out}")
    try:
        bp = BodyParts3D(args.bodyparts3d)
    except FileNotFoundError as exc:
        print(exc)
        return 1
    print(f"BodyParts3D: {len(bp.elements)} piezas.")

    print("Tablas de células y genes (Human Reference Atlas) y Uberon:")
    rows = []
    for organ in ASCTB_ORGANS:
        path = download(ASCTB_URL.format(organ=organ), CACHE / f"asct-b-vh-{organ}.csv")
        rows += read_asctb(path, organ)
    uberon = read_uberon_fma(download(UBERON_URL, CACHE / "uberon-edit.obo"))
    print(f"  {len(rows)} filas de {len(ASCTB_ORGANS)} tablas · {len(uberon)} equivalencias FMA–Uberon")

    def progress(key, data):
        if key == "anatomy.build.meshes":
            end = "\n" if data["done"] == data["total"] else ""
            print(f"\r  mallas: {data['done']}/{data['total']}", end=end, flush=True)

    report = build_anatomy(bp, rows, uberon, args.out, ratio=args.ratio, progress=progress)
    print(f"\nTerminado en {report.seconds:.0f} s: {report.pieces} piezas, "
          f"{report.triangles:,} triángulos, {report.bytes / 1e6:.1f} MB.".replace(",", "."))
    for key, info in report.systems.items():
        print(f"  {key:<13} {info['pieces']:>5} piezas  {info['bytes'] / 1e6:6.2f} MB")
    print(f"Células: {report.with_cells} piezas con las suyas (o las del órgano del que forman parte), "
          f"{report.with_tissue} con las del tejido general, {report.without} sin datos.")
    print(f"Datos: {ATTRIBUTION} · {ASCTB_CITATION} · Uberon (CC BY 3.0) · Cell Ontology (CC BY 4.0).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
