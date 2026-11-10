"""Prueba varias formas de pedir a ChEMBL las indicaciones de una enfermedad por su ID de MONDO.

    python tools/chembl_probe.py                  # anemia falciforme, MONDO:0011382
    python tools/chembl_probe.py MONDO:0007254    # otra enfermedad

Para cada forma enseña el código HTTP, cuánto tardó, cuántas indicaciones encontró y el
principio de la respuesta si es un error. Sin caché y sin tocar la base de datos.
"""

import json
import sys
import time
import urllib.error
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dnainvent.protein.clients import http_get  # noqa: E402

API = "https://www.ebi.ac.uk/chembl/api/data/drug_indication.json?"


def probe(label: str, url: str) -> None:
    started = time.monotonic()
    try:
        raw = http_get(url, timeout=90)        # como la importación
    except urllib.error.HTTPError as error:
        body = error.read()[:300]
        print(f"  {label:<34} HTTP {error.code} en {time.monotonic() - started:4.1f} s · {body!r}")
        return
    except Exception as error:          # tiempo agotado, sin red...
        print(f"  {label:<34} {type(error).__name__}: {error} ({time.monotonic() - started:.1f} s)")
        return
    elapsed = time.monotonic() - started
    try:
        data = json.loads(raw)
        total = (data.get("page_meta") or {}).get("total_count")
        first = (data.get("drug_indications") or [{}])[0]
        print(f"  {label:<34} OK en {elapsed:4.1f} s · total {total} · ej. {first.get('efo_id')} {first.get('efo_term')!r}")
    except ValueError:
        print(f"  {label:<34} OK en {elapsed:4.1f} s pero no es JSON: {raw[:200]!r}")


def main():
    mondo = sys.argv[1] if len(sys.argv) > 1 else "MONDO:0011382"
    number = mondo.split(":", 1)[-1]
    variants = [
        ("efo_id=MONDO:… (%3A)", API + urllib.parse.urlencode({"efo_id": mondo, "limit": 5})),
        ("efo_id=MONDO:… (: sin codificar)", API + f"efo_id={mondo}&limit=5"),
        ("efo_id__exact", API + urllib.parse.urlencode({"efo_id__exact": mondo, "limit": 5})),
        ("efo_id__iexact", API + urllib.parse.urlencode({"efo_id__iexact": mondo, "limit": 5})),
        ("efo_id__endswith=número", API + urllib.parse.urlencode({"efo_id__endswith": number, "limit": 5})),
        ("efo_id__contains=número", API + urllib.parse.urlencode({"efo_id__contains": number, "limit": 5})),
        ("efo_id=MONDO_… (guion bajo)", API + urllib.parse.urlencode({"efo_id": mondo.replace(':', '_'), "limit": 5})),
        ("mesh_id=D000755 (referencia)", API + urllib.parse.urlencode({"mesh_id": "D000755", "limit": 5})),
    ]
    print(f"Indicaciones de {mondo} en ChEMBL:")
    for label, url in variants:
        probe(label, url)
        time.sleep(0.5)


if __name__ == "__main__":
    main()
