"""Construye la base de datos de los genes de serie (el «catálogo» de la opción 9).

    python tools/build_catalog.py                       # los 14 genes -> ~/.dnainvent/catalog.db
    python tools/build_catalog.py --genes HBB,CFTR      # solo algunos
    python tools/build_catalog.py --max-variants 200    # prueba rápida (pocas variantes por gen)
    python tools/build_catalog.py --redo                # repetir también los pasos ya hechos
    python tools/build_catalog.py --out C:/ruta/catalog.db

Necesita conexión: NCBI Gene, ClinVar, UniProt, ChEMBL y Orphanet. Puede tardar bastante (los
genes como BRCA1 y BRCA2 tienen miles de variantes). Si se corta, vuelve a ejecutarlo: sigue
donde lo dejó.

Más rápido con una API key del NCBI (gratuita: cuenta del NCBI -> Settings -> API Key
Management), en la variable de entorno NCBI_API_KEY. Nunca la escribas en este archivo.

Si Orphanet no se puede descargar, bájate a mano en_product9_ages.xml y en_product9_prev.xml
de https://www.orphadata.com/ y pásale la carpeta con --orphadata-dir.
"""

import argparse
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dnainvent.database import open_database  # noqa: E402
from dnainvent.database.catalog import CATALOG_GENES, CatalogBuilder, StepResult  # noqa: E402
from dnainvent.database.orphanet import OrphanetProvider  # noqa: E402
from i18n import Translator  # noqa: E402

STEP_NAMES = {"gene": "gen (NCBI)", "variants": "variantes (ClinVar)",
              "sequences": "secuencias (RefSeq)", "gene_drugs": "fármacos de su proteína (ChEMBL)",
              "disease_drugs": "fármacos de sus enfermedades (ChEMBL)",
              "orphanet": "herencia y prevalencia (Orphanet)", "common_names": "nombres conocidos"}
MARKS = {"done": "✔", "skipped": "·", "failed": "✘"}

tr = Translator("es")


class Progress:
    """Como en la consola: los pasos por lotes se reescriben en la misma línea."""

    def __init__(self):
        self.same_line = None

    def __call__(self, key: str, data: dict):
        text = "    " + tr(key, **data)
        if "done" in data:
            if self.same_line == key:
                print("\r" + text, end="", flush=True)
            else:
                self.finish()
                print(text, end="", flush=True)
            self.same_line = key
        else:
            self.finish()
            print(text, flush=True)

    def finish(self):
        if self.same_line:
            print()
        self.same_line = None


def main() -> int:
    parser = argparse.ArgumentParser(description="Construye la base de datos de los genes de serie.")
    parser.add_argument("--out", default=str(Path.home() / ".dnainvent" / "catalog.db"),
                        help="archivo de la base de datos (por defecto ~/.dnainvent/catalog.db)")
    parser.add_argument("--genes", help="genes separados por comas (por defecto, los 14 de serie)")
    parser.add_argument("--max-variants", type=int, help="como mucho tantas variantes por gen (pruebas)")
    parser.add_argument("--redo", action="store_true", help="repetir también los pasos ya hechos")
    parser.add_argument("--orphadata-dir", help="carpeta con los XML de Orphadata ya descargados")
    args = parser.parse_args()

    genes = [g.strip().upper() for g in args.genes.split(",") if g.strip()] if args.genes else CATALOG_GENES
    out = Path(args.out).expanduser()
    print(f"Catálogo de DNAinvent: {len(genes)} genes -> {out}")
    if os.environ.get("NCBI_API_KEY"):
        print("NCBI_API_KEY encontrada: 10 consultas por segundo.")
    else:
        print("Sin NCBI_API_KEY: 3 consultas por segundo (funciona igual, pero más despacio).")
    if args.max_variants:
        print(f"Modo de prueba: como mucho {args.max_variants} variantes por gen.")
    progress = Progress()

    def on_gene(symbol: str):
        progress.finish()
        print("\nTodo el catálogo" if symbol == "*" else f"\n{symbol}", flush=True)

    def announce(result: StepResult):
        progress.finish()
        line = f"  {MARKS[result.status]} {STEP_NAMES[result.step]}"
        if result.status == "skipped":
            line += " (ya estaba)"
        if result.detail:
            line += f": {result.detail}"
        print(line, flush=True)

    started = time.monotonic()
    conn = open_database(out)
    try:
        builder = CatalogBuilder(conn, progress=progress, announce=announce, on_gene=on_gene,
                                 max_variants=args.max_variants,
                                 orphanet=OrphanetProvider(local_dir=args.orphadata_dir))
        report = builder.run(genes, redo=args.redo)
    except KeyboardInterrupt:
        progress.finish()
        print("\nInterrumpido. Vuelve a ejecutarlo para seguir donde lo has dejado.")
        return 130
    finally:
        conn.close()

    minutes = (time.monotonic() - started) / 60
    if report.removed_genes:
        print("\nGenes borrados porque no son del catálogo: " + ", ".join(report.removed_genes))
    size = out.stat().st_size / 1_000_000
    print(f"\nTerminado en {minutes:.1f} min. Tamaño: {size:.1f} MB.")
    if report.missing_names:
        print("Nombres conocidos que no se han encontrado (revísalos): "
              + ", ".join(f"{name} ({gene})" for gene, name in report.missing_names))
    if report.failed:
        print(f"\nPasos fallidos: {len(report.failed)}. Vuelve a ejecutarlo más tarde para repetirlos:")
        for step in report.failed:
            print(f"  ✘ {step.symbol} · {STEP_NAMES[step.step]}: {step.detail}")
        return 1
    print("Todos los pasos están hechos.")
    print("Datos: NCBI Gene, RefSeq y ClinVar (dominio público) · ChEMBL (CC BY-SA 3.0) · "
          "UniProt (CC BY 4.0) · Orphanet/Orphadata (CC BY 4.0).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
