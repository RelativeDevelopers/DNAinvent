"""¿Qué es la proteína del usuario? Los cuatro casos, a partir de los números de BLAST.

BLAST compara la proteína del usuario (la "consulta") con una base de proteínas reales y
devuelve sus parecidos ("hits"). De cada hit importan tres números:

    identidad    % de aminoácidos iguales en el tramo alineado
    cobertura    qué parte de la consulta (y de la proteína real) entra en ese tramo
    E-value      cuántos parecidos así de buenos saldrían por pura casualidad
                 (1e-50: imposible que sea casualidad; 1: lo normal por azar)

La identidad sola engaña: 12 aminoácidos idénticos de una proteína de 800 no significan nada.
Por eso un hit solo cuenta si su E-value es bajo y cubre buena parte de la consulta.

    Caso               Condición (valores por defecto, en Thresholds)
    identical          100 % de identidad, cubriendo entera la consulta y la proteína real
    very_similar       identidad >= 90 %: la misma proteína con mutaciones
    distant_relative   identidad >= 30 %: comparten antepasado (p. ej. otra especie)
    invented           ningún hit válido

Por debajo de ~30 % está la "zona crepuscular" (twilight zone, Rost 1999): ahí ya no se
distingue un parentesco real de una coincidencia, así que se trata como inventada.

Este módulo no se conecta a nada: recibe los números ya calculados (el cliente de BLAST
los sacará de la respuesta del EBI) y así se puede probar sin red.
"""

from dataclasses import dataclass
from enum import Enum


class Category(str, Enum):
    IDENTICAL = "identical"
    VERY_SIMILAR = "very_similar"
    DISTANT_RELATIVE = "distant_relative"
    INVENTED = "invented"

    @property
    def i18n_key(self) -> str:
        """Clave base de los textos: protein.categories.<caso>.name / .description"""
        return f"protein.categories.{self.value}"


# De mejor a peor: si hay varios hits válidos, manda el de mejor caso.
_RANK = {Category.IDENTICAL: 0, Category.VERY_SIMILAR: 1, Category.DISTANT_RELATIVE: 2}
_FULL = 100.0 - 1e-9       # para comparar porcentajes calculados con decimales


@dataclass(frozen=True)
class Hit:
    """Un parecido encontrado por BLAST (solo los números que hacen falta para decidir).

    accession         identificador de UniProt, p. ej. "P69905" (hemoglobina α humana)
    identity          % de aminoácidos idénticos en el tramo alineado (0-100)
    query_coverage    % de la proteína del usuario que entra en el tramo (0-100)
    subject_coverage  % de la proteína real que entra en el tramo (0-100)
    evalue            E-value (>= 0)
    preferred         desempate: si dos hits empatan en todo, gana el preferido (la
                      proteína humana: en el caso de la hemoglobina, la del chimpancé y la
                      humana son idénticas y BLAST puede poner primero cualquiera)
    """

    accession: str
    identity: float
    query_coverage: float
    subject_coverage: float
    evalue: float
    preferred: bool = False

    def __post_init__(self):
        for name in ("identity", "query_coverage", "subject_coverage"):
            value = getattr(self, name)
            if not 0 <= value <= 100:
                raise ValueError(f"{name} debe estar entre 0 y 100: {value}")
        if self.evalue < 0:
            raise ValueError(f"evalue no puede ser negativo: {self.evalue}")


@dataclass(frozen=True)
class Thresholds:
    """Umbrales de la clasificación. Cambiarlos aquí cambia toda la aplicación."""

    very_similar_identity: float = 90.0
    distant_identity: float = 30.0
    min_query_coverage: float = 70.0
    max_evalue: float = 1e-5

    def __post_init__(self):
        if not 0 < self.distant_identity < self.very_similar_identity < 100:
            raise ValueError("Debe cumplirse 0 < distant_identity < very_similar_identity < 100")


@dataclass(frozen=True)
class Classification:
    category: Category
    hit: Hit | None             # el hit que decide el caso (None si es inventada)

    @property
    def is_known(self) -> bool:
        """True si hay una proteína real con la que comparar."""
        return self.hit is not None


def category_of(hit: Hit, thresholds: Thresholds = Thresholds()) -> Category | None:
    """Caso al que llevaría este hit por sí solo, o None si no es fiable."""
    if hit.evalue > thresholds.max_evalue or hit.query_coverage < thresholds.min_query_coverage:
        return None
    if hit.identity >= _FULL and hit.query_coverage >= _FULL and hit.subject_coverage >= _FULL:
        return Category.IDENTICAL
    if hit.identity >= thresholds.very_similar_identity:
        return Category.VERY_SIMILAR
    if hit.identity >= thresholds.distant_identity:
        return Category.DISTANT_RELATIVE
    return None


def classify(hits: list[Hit], thresholds: Thresholds = Thresholds()) -> Classification:
    """Elige el caso a partir de todos los hits de una búsqueda.

    Gana el hit con mejor caso; a igualdad, el de más identidad, después el de menor
    E-value, después el de más cobertura de la consulta y después el preferido. Si aún
    empatan, el primero de la lista (el orden de BLAST).
    """
    candidates = [(category, hit) for hit in hits
                  if (category := category_of(hit, thresholds)) is not None]
    if not candidates:
        return Classification(Category.INVENTED, None)
    category, hit = min(candidates, key=lambda c: (_RANK[c[0]], -c[1].identity, c[1].evalue,
                                                   -c[1].query_coverage, not c[1].preferred))
    return Classification(category, hit)
