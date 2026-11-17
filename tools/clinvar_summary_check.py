"""Diagnóstico: qué dice ClinVar de una variante y qué hace DNAinvent con ella.

    python tools/clinvar_summary_check.py rs5030868 G6PD     # por rs (o cualquier búsqueda)
    python tools/clinvar_summary_check.py 100057 G6PD        # por ID de variante de ClinVar

Para cada registro: su título, su tipo, las variantes que contiene, sus genes, si la
importación lo guardaría (y por qué no) y si sale al buscar las variantes del gen, que es lo
que hace la importación. Sirve para entender por qué una variante conocida no está en el
catálogo.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dnainvent.database.clinvar import ClinVarProvider, is_combined, parse_summary  # noqa: E402


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return
    query = sys.argv[1]
    symbol = sys.argv[2].upper() if len(sys.argv) > 2 else None
    client = ClinVarProvider()
    if query.isdigit():
        ids = [query]
    else:
        ids = client.json("esearch.fcgi", {"db": "clinvar", "term": query, "retmax": 20}) \
            .get("esearchresult", {}).get("idlist", [])
        print(f"Búsqueda «{query}»: {len(ids)} registros {ids}")
    if not ids:
        return
    in_gene = set(client.search_variant_ids(symbol)) if symbol else set()
    if symbol:
        print(f"Búsqueda «{symbol}[gene]» (la de la importación): {len(in_gene)} registros")

    for uid, doc in client.fetch_summaries(ids):
        print(f"\n=== ID {uid} · {doc.get('accession')} ===")
        print(f"  título:    {doc.get('title')}")
        print(f"  tipo:      {doc.get('obj_type')}")
        for item in doc.get("variation_set") or []:
            print(f"  variante:  {item.get('variation_name')}  [{item.get('variant_type')}]")
        print(f"  genes:     {[g.get('symbol') for g in doc.get('genes') or []]}")
        germline = doc.get("germline_classification") or {}
        print(f"  clasific.: {germline.get('description')} ({germline.get('review_status')})")
        parsed = parse_summary(uid, doc)
        if parsed is None:
            why = "es una combinación de variantes" if is_combined(doc) else "no está descrita sobre un ARNm"
            print(f"  DNAinvent: NO la guardaría ({why})")
        else:
            print(f"  DNAinvent: {parsed['accession']}.{parsed['version']} gen {parsed['gene']} "
                  f"{parsed['hgvs_c']} {parsed['hgvs_p']} tipo {parsed['variant_type']}")
            if symbol and parsed["gene"].upper() != symbol:
                print(f"             NO la guardaría en {symbol}: el nombre dice {parsed['gene']}")
        if symbol:
            print(f"  ¿sale en «{symbol}[gene]»?  {'sí' if uid in in_gene else 'NO'}")


if __name__ == "__main__":
    main()
