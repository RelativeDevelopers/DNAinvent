"""Pequeñas utilidades compartidas por los módulos de la base de datos."""

from datetime import datetime, timezone


def now_iso() -> str:
    """Fecha y hora actuales en UTC, en ISO 8601 (así se guardan en la base de datos)."""
    return datetime.now(timezone.utc).isoformat(timespec="seconds")
