"""Lector y escritor de archivos PDB (solo lo necesario, sin dependencias).

El formato PDB (1976) es de columnas fijas: cada átomo es una línea de 80 caracteres y cada
dato está siempre en las mismas columnas. Por eso se lee cortando la línea, sin buscar
espacios (que a veces faltan, cuando dos números largos quedan pegados):

    ATOM      2  CA  MET A   1      51.653 -89.304   8.833  1.00 10.38           C
    |    |      |    |   | |        |       |       |       |     |               |
    tipo número nombre res cadena   x       y       z      ocup. B (*)        elemento

(*) En las predicciones (AlphaFold, ESMFold) la columna B no es la "agitación térmica" de
un experimento sino la confianza de cada átomo, el pLDDT (0-100).

Solo se leen las líneas ATOM y HETATM del primer modelo (las estructuras de RMN traen
varios). Si un átomo tiene posiciones alternativas (altLoc), se queda la primera.

AlphaFold DB y ESMFold devuelven PDB, así que con esto basta. Las estructuras experimentales
muy grandes del PDB solo existen en mmCIF; se añadirá un lector cuando haga falta.
"""

from dataclasses import dataclass, field
from pathlib import Path

from ..constants import AMINO_ACID_CODES

# "ALA" -> "A". Lo que no está aquí (ligandos, agua, aminoácidos raros) es "X".
THREE_TO_ONE = {three.upper(): one for one, three in AMINO_ACID_CODES.items()}


def normalize_name(name: str) -> str:
    """Nombres antiguos de hidrógenos (PDB v2, los que usa PyMOL) -> actuales (v3).

        "2HB" -> "HB2"     "1HD2" -> "HD21"     "HA" -> "HA"
    """
    return name[1:] + name[0] if name[:1].isdigit() else name


@dataclass(frozen=True)
class Atom:
    name: str              # "CA", "HB2"...
    resn: str              # residuo en 3 letras: "MET"
    chain: str
    resi: int              # número de residuo tal como viene en el archivo
    x: float               # Å
    y: float
    z: float
    element: str           # "C", "N", "S"...
    bfactor: float = 0.0   # en predicciones: pLDDT
    icode: str = ""        # código de inserción (raro; numeraciones como 52, 52A, 53)
    hetero: bool = False   # True en líneas HETATM (ligandos, agua...)
    occupancy: float = 1.0

    @property
    def xyz(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.z)


@dataclass
class Residue:
    chain: str
    resi: int
    icode: str
    resn: str
    atoms: dict[str, Atom] = field(default_factory=dict)    # nombre -> átomo

    @property
    def letter(self) -> str:
        return THREE_TO_ONE.get(self.resn, "X")


@dataclass
class Structure:
    atoms: list[Atom]
    title: str = ""

    @property
    def is_prediction(self) -> bool:
        """True si viene de AlphaFold o ESMFold (la columna B es entonces el pLDDT)."""
        title = self.title.upper()
        return "ALPHAFOLD" in title or "ESMFOLD" in title

    def chains(self) -> list[str]:
        """Cadenas con aminoácidos, en el orden del archivo."""
        seen = []
        for atom in self.atoms:
            if not atom.hetero and atom.chain not in seen:
                seen.append(atom.chain)
        return seen

    def residues(self, chain: str | None = None) -> list[Residue]:
        """Aminoácidos (líneas ATOM) de una cadena, en orden. Por defecto, la primera."""
        chains = self.chains()
        if not chains:
            return []
        chain = chains[0] if chain is None else chain
        out, index = [], {}
        for atom in self.atoms:
            if atom.hetero or atom.chain != chain:
                continue
            key = (atom.resi, atom.icode)
            if key not in index:
                index[key] = len(out)
                out.append(Residue(chain, atom.resi, atom.icode, atom.resn))
            out[index[key]].atoms.setdefault(atom.name, atom)
        return out

    def sequence(self, chain: str | None = None) -> str:
        return "".join(residue.letter for residue in self.residues(chain))


def parse_pdb(text: str, model: int = 1) -> Structure:
    """Texto de un archivo PDB -> Structure (solo el modelo indicado, por defecto el 1)."""
    atoms, title = [], []
    current_model = 1
    alt_chosen = {}                       # (cadena, residuo, icode, nombre) -> altLoc elegido
    for line in text.splitlines():
        record = line[:6]
        if record == "TITLE ":
            title.append(line[10:80].strip())
        elif record == "MODEL ":
            current_model = int(line[10:14])
        elif record == "ENDMDL" and current_model >= model:
            break
        elif record in ("ATOM  ", "HETATM") and current_model == model:
            atom = _parse_atom(line)
            alt = line[16:17]
            key = (atom.chain, atom.resi, atom.icode, atom.name)
            if alt.strip():
                if alt_chosen.setdefault(key, alt) != alt:
                    continue              # otra posición alternativa del mismo átomo
            atoms.append(atom)
    return Structure(atoms, " ".join(title))


def read_pdb(path: str | Path, model: int = 1) -> Structure:
    return parse_pdb(Path(path).read_text(encoding="utf-8", errors="replace"), model)


def _parse_atom(line: str) -> Atom:
    line = line.ljust(80)
    name = normalize_name(line[12:16].strip())
    element = line[76:78].strip()
    if not element:                       # archivos antiguos sin columna de elemento
        element = next(c for c in name if c.isalpha())
    return Atom(
        name=name,
        resn=line[17:20].strip(),
        chain=line[21:22].strip() or "A",
        resi=int(line[22:26]),
        icode=line[26:27].strip(),
        x=float(line[30:38]), y=float(line[38:46]), z=float(line[46:54]),
        occupancy=float(line[54:60] or 1.0),
        bfactor=float(line[60:66] or 0.0),
        element=element.capitalize(),
        hetero=line[:6] == "HETATM",
    )


def write_pdb(atoms: list[Atom], title: str = "") -> str:
    """Átomos -> texto PDB (para abrirlo en PyMOL, ChimeraX, Mol*...)."""
    lines = [f"TITLE     {title[:70]}"] if title else []
    for serial, a in enumerate(atoms, start=1):
        # Convención de columnas: los elementos de 1 letra empiezan en la columna 14
        name = f" {a.name:<3}" if len(a.element) == 1 and len(a.name) < 4 else f"{a.name:<4}"
        record = "HETATM" if a.hetero else "ATOM  "
        lines.append(
            f"{record}{serial % 100000:5d} {name} {a.resn:>3} {a.chain[:1]:1}{a.resi:4d}"
            f"{a.icode[:1]:1}   {a.x:8.3f}{a.y:8.3f}{a.z:8.3f}{a.occupancy:6.2f}"
            f"{a.bfactor:6.2f}          {a.element:>2}"
        )
    lines.append("END")
    return "\n".join(lines) + "\n"
