"""Los elementos del ADN (H, C, N, O y P) y de las proteínas (H, C, N, O y S), con los datos
de AtomicWorld.

De la tabla periódica de AtomicWorld salen el radio de van der Waals, el color convencional
(Jmol: no son colores reales, los átomos son mucho más pequeños que la longitud de onda de
la luz), la configuración electrónica y el isótopo más abundante (H-1, C-12, N-14, O-16, P-31).
"""
from dataclasses import dataclass
from functools import lru_cache

from atomicworld import atoms
from atomicworld.repository import IsotopeTable, PeriodicTable

DNA_ELEMENTS = ("H", "C", "N", "O", "P")
PROTEIN_ELEMENTS = ("H", "C", "N", "O", "S")     # S: cisteína y metionina


@dataclass(frozen=True)
class AtomType:
    symbol: str
    element: object          # atomicworld.models.Element
    isotope: object          # atomicworld.models.Isotope (el más abundante)
    shells: tuple            # subniveles ocupados (atomicworld.atoms.Subshell)

    @property
    def number(self) -> int:
        return self.element.number

    @property
    def electrons(self) -> int:
        return self.element.electrons

    @property
    def vdw_nm(self) -> float:
        return self.element.vdw_radius_pm / 1000

    @property
    def color(self) -> tuple[float, float, float]:
        """'#3050f8' -> (0.19, 0.31, 0.97)"""
        h = self.element.color.lstrip("#")
        return tuple(int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))


@lru_cache(maxsize=1)
def _tables():
    return PeriodicTable(), IsotopeTable()


@lru_cache(maxsize=None)
def atom_type(symbol: str) -> AtomType:
    table, isotopes = _tables()
    element = table.find(symbol)
    return AtomType(symbol, element, isotopes.default(element), atoms.subshells(element))
