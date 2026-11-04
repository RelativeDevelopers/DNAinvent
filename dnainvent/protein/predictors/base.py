"""Interfaz común de los predictores de estructura.

    class MiPredictor:
        name = "Mi predictor"        # para mostrarlo y citarlo
        max_length = 400            # aminoácidos como máximo
        def predict(self, protein: str) -> Structure: ...

Reglas:
  - Devuelve una Structure con el pLDDT (0-100) en la columna B de cada átomo.
  - Si la proteína es demasiado larga: InvalidParameterError("prediction_too_long").
  - Si el servicio no responde: ServiceError("service_unavailable", service=name).
  - Si responde pero no puede con esta proteína: ServiceError("prediction_failed", ...).
"""
from typing import Protocol

from ..pdb import Structure


class Predictor(Protocol):
    name: str
    max_length: int

    def predict(self, protein: str) -> Structure: ...
