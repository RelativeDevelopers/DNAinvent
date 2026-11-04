"""DNAinvent: herramientas para trabajar con secuencias de ADN y ARN.

Este paquete solo contiene lógica: no pregunta nada ni imprime nada,
así cualquier interfaz (consola, ventana, web) puede usarlo.
"""

from .errors import DNAinventError, InvalidParameterError, InvalidSequenceError, ServiceError

__version__ = "2.1.0-dev"
__all__ = ["DNAinventError", "InvalidParameterError", "InvalidSequenceError", "ServiceError"]
