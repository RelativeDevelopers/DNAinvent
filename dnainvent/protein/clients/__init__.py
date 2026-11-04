"""Clientes de servicios externos (AlphaFold DB, BLAST del EBI; más adelante ESMFold).

Todos usan urllib (biblioteca estándar) y reciben la función de descarga como parámetro,
para poder probarlos sin red. Los fallos se convierten en errores con código
(protein_not_found, service_unavailable...) que cada interfaz traduce.
"""
import urllib.parse
import urllib.request

from ... import __version__

USER_AGENT = f"DNAinvent/{__version__} (+https://github.com/aluskabarbera/DNAinvent)"
TIMEOUT = 30          # segundos


def http_get(url: str, timeout: float = TIMEOUT) -> bytes:
    """Descarga una URL. Lanza urllib.error.HTTPError (404...) o URLError/OSError (sin red)."""
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


def http_post(url: str, data: dict, timeout: float = TIMEOUT) -> bytes:
    """Envía un formulario (application/x-www-form-urlencoded) y devuelve la respuesta."""
    body = urllib.parse.urlencode(data).encode("utf-8")
    request = urllib.request.Request(url, data=body, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


def http_post_text(url: str, text: str, timeout: float = TIMEOUT) -> bytes:
    """Envía un texto tal cual como cuerpo de la petición (como `curl --data`)."""
    request = urllib.request.Request(url, data=text.encode("utf-8"), headers={
        "User-Agent": USER_AGENT, "Content-Type": "application/x-www-form-urlencoded"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()
