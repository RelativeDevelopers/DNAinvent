"""
Prueba de la fase 2: descarga (o lee) una estructura, le añade los hidrógenos y la resume.

    python tools/protein_info.py P68871                       # de AlphaFold DB (con caché)
    python tools/protein_info.py --file proteina.pdb          # un archivo local
    python tools/protein_info.py P68871 --save hbb_con_h.pdb  # guarda la estructura con H

El archivo de --save se puede abrir en PyMOL, ChimeraX o https://molstar.org/viewer/
para ver los hidrógenos calculados.

Herramienta para desarrollar: los textos van en español, sin traducciones.
"""
import argparse
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dnainvent import DNAinventError                                   # noqa: E402
from dnainvent.protein import confidence                                # noqa: E402
from dnainvent.protein.clients.afdb import AlphaFoldDB                  # noqa: E402
from dnainvent.protein.hydrogens import add_hydrogens                   # noqa: E402
from dnainvent.protein.pdb import read_pdb, write_pdb                   # noqa: E402
from dnainvent.scene.protein import build_protein                       # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[1].strip())
    parser.add_argument("accession", nargs="?", help="accesión de UniProt, p. ej. P68871")
    parser.add_argument("--file", help="leer un archivo PDB local en vez de descargar")
    parser.add_argument("--save", help="guardar la estructura con hidrógenos en este PDB")
    args = parser.parse_args()
    if not args.accession and not args.file:
        parser.error("indica una accesión o --file")

    try:
        if args.file:
            structure = read_pdb(args.file)
            print(f"Archivo: {args.file}")
        else:
            db = AlphaFoldDB()
            entry, structure = db.structure(args.accession)
            print(f"Modelo:   {entry.model_id} (versión {entry.version})")
            print(f"Caché:    {db.cache.root}")
            if entry.sequence and entry.sequence != structure.sequence():
                print("AVISO: la secuencia del PDB no coincide con la de la ficha")
    except DNAinventError as error:
        sys.exit(f"Error: {error}")

    residues = structure.residues()
    sequence = structure.sequence()
    print(f"Título:   {structure.title}")
    print(f"Cadenas:  {', '.join(structure.chains())}  (se usa la primera)")
    print(f"Longitud: {len(sequence)} aminoácidos")
    print(f"Secuencia: {sequence[:60]}{'…' if len(sequence) > 60 else ''}")

    if structure.is_prediction:
        plddt = [r.atoms["CA"].bfactor for r in residues if "CA" in r.atoms]
        print(f"pLDDT medio: {sum(plddt) / len(plddt):.1f}")
        bands = Counter(confidence.band(p) for p in plddt)
        for _, key, _ in confidence.BANDS:
            print(f"   {key:<10} {bands[key]:4d} aminoácidos")

    result = add_hydrogens(residues)
    heavy = sum(1 for a in structure.atoms if not a.hetero and a.element != "H")
    print(f"Átomos pesados: {heavy}   H calculados: {result.added}")
    if result.skipped:
        print(f"H sin calcular (faltan átomos de referencia): {result.skipped[:10]}")

    try:
        molecule = build_protein(structure)
    except ImportError as error:              # build_protein necesita NumPy
        print(f"(sin NumPy no se construye la molécula 3D: {error})")
    else:
        elements = Counter(molecule.elements)
        print(f"Molécula 3D: {molecule.n_atoms} átomos {dict(elements)}, "
              f"{len(molecule.bonds)} enlaces")

    if args.save:
        atoms = [atom for residue in result.residues for atom in residue.atoms.values()]
        Path(args.save).write_text(write_pdb(atoms, structure.title + " + H (DNAinvent)"),
                                   encoding="utf-8")
        print(f"Guardado: {args.save}")


if __name__ == "__main__":
    main()
