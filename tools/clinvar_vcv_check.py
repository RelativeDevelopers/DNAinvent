"""Diagnóstico: cómo responde ClinVar a la petición de clasificaciones por enfermedad (XML VCV).

    python tools/clinvar_vcv_check.py            # TP53
    python tools/clinvar_vcv_check.py HBB

Prueba la misma petición que hace la importación con 1, 3 y 25 variantes del gen y dice, para
cada una, si llegó, cuánto ocupa, cuántos registros trae y cuánto tardó. Guarda las respuestas
en clinvar_vcv_*.xml (en la carpeta desde la que lo ejecutes).
"""

import sys
import time
import urllib.error
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dnainvent.database.clinvar import ClinVarProvider, parse_vcv  # noqa: E402


def main():
    symbol = sys.argv[1] if len(sys.argv) > 1 else "TP53"
    client = ClinVarProvider()
    ids = client.search_variant_ids(symbol, max_variants=25)
    print(f"{symbol}: primeros IDs {ids[:5]}")
    for n in (1, 3, 25):
        chunk = ids[:n]
        url = client.url("efetch.fcgi", {"db": "clinvar", "rettype": "vcv", "is_variationid": "true",
                                         "id": ",".join(chunk)})
        print(f"\n=== {n} variante(s) ===")
        started = time.monotonic()
        try:
            raw = client.http(url)
        except urllib.error.HTTPError as error:
            print(f"  HTTP {error.code}: {error.read()[:300]!r}")
            continue
        except Exception as error:      # tiempo agotado, sin red...
            print(f"  ERROR {type(error).__name__}: {error}  ({time.monotonic() - started:.1f} s)")
            continue
        print(f"  {len(raw):,} bytes en {time.monotonic() - started:.1f} s")
        path = Path(f"clinvar_vcv_{n}.xml")
        path.write_bytes(raw)
        try:
            root = ET.fromstring(raw)
        except ET.ParseError as error:
            print(f"  no es XML válido: {error}; empieza por {raw[:200]!r}")
            continue
        archives = list(root.iter("VariationArchive"))
        print(f"  raíz <{root.tag}>, hijos {[child.tag for child in root][:5]}")
        print(f"  VariationArchive: {len(archives)}  RCVAccession: {len(list(root.iter('RCVAccession')))}"
              f"  ClassifiedCondition: {len(list(root.iter('ClassifiedCondition')))}")
        for archive in archives[:3]:
            print(f"    atributos: {dict(list(archive.attrib.items())[:6])}")
        returned = {int(a.get("VariationID", 0)) for a in archives}
        asked = {int(i) for i in chunk}
        print(f"  pedidas que han llegado: {len(returned & asked)} de {len(asked)}; "
              f"no pedidas: {len(returned - asked)}")
        conditions = {v: c for v, c in parse_vcv(raw).items() if v in asked}
        print(f"  DNAinvent lee: {sum(len(v) for v in conditions.values())} condiciones en "
              f"{sum(1 for v in conditions.values() if v)} de {len(chunk)} variantes")
        print(f"  guardado en {path.resolve()}")


if __name__ == "__main__":
    main()
