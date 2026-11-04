"""Comparar la proteína del usuario con la natural: secuencia y forma.

1. Qué aminoácido corresponde a cuál. Se usa el alineamiento que ya hizo BLAST (lo que el
   usuario ve en pantalla); si no lo hay, se alinea aquí (alignment.py, BLOSUM62).
2. Diferencias de secuencia, con el nombre habitual de una mutación:
       Glu7Val = en la posición 7 de la natural hay Glu (glutamato) y en la tuya, Val (valina)
   La numeración es la de UniProt (empieza en la metionina inicial). Por eso la mutación de
   la anemia falciforme aparece como Glu7Val, aunque en medicina se llame Glu6Val (la
   numeración clásica no cuenta esa metionina).
3. Diferencias de forma: se superponen los carbonos α (Cα, uno por aminoácido) de las
   parejas y se mide el RMSD y la distancia de cada pareja después de superponer.
"""

from dataclasses import dataclass, field, replace

from ..alignment import GAP, align
from ..constants import AMINO_ACID_CODES
from .clients.blast import BlastHit
from .pdb import Structure
from .superpose import Superposition, superpose


@dataclass(frozen=True)
class Difference:
    natural_position: int     # desde 1, numeración de la proteína natural
    user_position: int        # desde 1, numeración de la proteína del usuario
    natural: str              # aminoácido (1 letra)
    user: str

    @property
    def label(self) -> str:
        """"Glu7Val" """
        return (f"{AMINO_ACID_CODES.get(self.natural, 'Xaa')}{self.natural_position}"
                f"{AMINO_ACID_CODES.get(self.user, 'Xaa')}")


@dataclass
class Comparison:
    pairs: list[tuple[int, int]]                   # (usuario, natural), desde 0
    differences: list[Difference]
    insertions: int                                # aminoácidos de más en la del usuario
    deletions: int                                 # aminoácidos que le faltan
    superposition: Superposition | None = None     # predicción -> natural
    deviations: dict[int, float] = field(default_factory=dict)   # natural (desde 0) -> Å

    @property
    def identity(self) -> float:
        """% de parejas con el mismo aminoácido."""
        return (len(self.pairs) - len(self.differences)) / len(self.pairs) * 100

    @property
    def rmsd(self) -> float | None:
        return self.superposition.rmsd if self.superposition else None

    @property
    def conserved(self) -> list[int]:
        """Posiciones de la natural (desde 0) que no han cambiado."""
        changed = {d.natural_position - 1 for d in self.differences}
        return [j for _, j in self.pairs if j not in changed]


def pairs_from_blast(hit: BlastHit) -> list[tuple[int, int]]:
    """Parejas (usuario, natural) desde 0, leyendo el alineamiento de BLAST columna a columna."""
    i, j, out = hit.query_from - 1, hit.hit_from - 1, []
    for a, b in zip(hit.query_aligned, hit.hit_aligned):
        if a != GAP and b != GAP:
            out.append((i, j))
        i += a != GAP
        j += b != GAP
    return out


def pairs_by_alignment(user: str, natural: str) -> list[tuple[int, int]]:
    return align(user, natural, kind="protein").pairs()


def compare(user: str, natural: str, hit: BlastHit | None = None,
            natural_structure: Structure | None = None,
            predicted: Structure | None = None) -> Comparison:
    """Compara la proteína del usuario con la natural.

    hit: el parecido de BLAST (si se pasa y trae el alineamiento, se usa ese).
    natural_structure / predicted: si se pasan las dos, también se compara la forma.
    """
    usable = hit is not None and hit.query_aligned and hit.hit_aligned
    pairs = pairs_from_blast(hit) if usable else pairs_by_alignment(user, natural)
    pairs = [(i, j) for i, j in pairs if i < len(user) and j < len(natural)]
    differences = [Difference(j + 1, i + 1, natural[j], user[i])
                   for i, j in pairs if user[i] != natural[j]]
    insertions = deletions = 0
    if usable:
        insertions = hit.hit_aligned.count(GAP)
        deletions = hit.query_aligned.count(GAP)
    else:
        result = align(user, natural, kind="protein")
        insertions, deletions = result.second.count(GAP), result.first.count(GAP)
    comparison = Comparison(pairs, differences, insertions, deletions)

    if natural_structure is not None and predicted is not None:
        _compare_shapes(comparison, natural_structure, predicted)
    return comparison


def _compare_shapes(comparison: Comparison, natural: Structure, predicted: Structure) -> None:
    natural_ca = [r.atoms.get("CA") for r in natural.residues()]
    predicted_ca = [r.atoms.get("CA") for r in predicted.residues()]
    both = [(i, j) for i, j in comparison.pairs
            if i < len(predicted_ca) and j < len(natural_ca)
            and predicted_ca[i] is not None and natural_ca[j] is not None]
    if len(both) < 3:
        return
    mobile = [predicted_ca[i].xyz for i, _ in both]
    reference = [natural_ca[j].xyz for _, j in both]
    fit = superpose(mobile, reference)
    comparison.superposition = fit
    comparison.deviations = {j: _distance(fit.apply(m), r)
                             for (_, j), m, r in zip(both, mobile, reference)}


def _distance(a, b) -> float:
    return sum((x - y) ** 2 for x, y in zip(a, b)) ** 0.5


def moved(structure: Structure, fit: Superposition) -> Structure:
    """Copia de la estructura con todos sus átomos movidos (para guardarla superpuesta)."""
    atoms = []
    for atom in structure.atoms:
        x, y, z = fit.apply(atom.xyz)
        atoms.append(replace(atom, x=x, y=y, z=z))
    return Structure(atoms, structure.title)
