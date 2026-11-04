"""Formateo de longitudes y tiempos para el HUD, con el separador decimal del idioma."""

_LENGTH_UNITS = [(1e-3, "mm"), (1e-6, "µm"), (1e-9, "nm"), (1e-12, "pm"),
                 (1e-15, "fm"), (1e-18, "am")]
_TIME_UNITS = [(1, "s"), (1e-3, "ms"), (1e-6, "µs"), (1e-9, "ns"), (1e-12, "ps"),
               (1e-15, "fs"), (1e-18, "as"), (1e-21, "zs"), (1e-24, "ys")]


def _fmt(value: float, units, dec: str) -> str:
    if abs(value) < 1e-30:
        return "0 m" if units is _LENGTH_UNITS else "0 s"     # "0 m/s", no "0/s"
    for factor, name in units:
        if abs(value) >= factor * 0.999:
            v = value / factor
            txt = f"{v:.0f}" if v >= 100 else f"{v:.1f}" if v >= 10 else f"{v:.2f}"
            return txt.replace(".", dec) + " " + name
    factor, name = units[-1]
    return f"{value / factor:.3f}".replace(".", dec) + " " + name


def length(meters: float, dec: str = ".") -> str:
    """1.2e-10 -> '120 pm'"""
    return _fmt(meters, _LENGTH_UNITS, dec)


def duration(seconds: float, dec: str = ".") -> str:
    """Tiempos cortos con prefijos del SI: 3.5e-22 -> '0.350 zs'"""
    return _fmt(seconds, _TIME_UNITS, dec)
