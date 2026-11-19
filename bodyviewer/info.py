"""Lo que el visor enseña de la pieza elegida (un resumen de la ficha de la consola)."""

CELLS_SHOWN = 5
GENES_SHOWN = 10


def shown_name(item: dict) -> str:
    return item.get("label") or item["name"]


def selection_info(card: dict | None, tr) -> dict | None:
    """{'title', 'subtitle', 'lines': [(texto, 'normal' | 'dim' | 'gold')]} para el HUD."""
    if card is None:
        return None
    lines = []
    facts = [tr(f"anatomy.systems.{card['system']}")] if card.get("system") else []
    if card.get("size_cm"):
        facts.append(tr("anatomy.card.size_value", cm=tr.decimal(card["size_cm"], 1, trim=True)))
    if facts:
        lines.append((" · ".join(facts), "dim"))
    if card["part_of"]:
        lines.append((f"{tr('anatomy.card.part_of')}: " + ", ".join(shown_name(p) for p in card["part_of"][:2]), "dim"))
    cells = card["cells"]
    if not cells["cells"]:
        lines.append((tr("anatomy.viewer.no_cells"), "normal"))
    else:
        names = [shown_name(c) for c in cells["cells"][:CELLS_SHOWN]]
        extra = len(cells["cells"]) - CELLS_SHOWN
        text = tr("anatomy.viewer.cells", list=", ".join(names) + (f" (+{extra})" if extra > 0 else ""))
        if cells["source"] == "tissue":
            text += " " + tr("anatomy.viewer.tissue", tissue=tr(f"anatomy.tissues.{cells['tissue']}"))
        lines.append((text, "normal"))
        genes = [("★" if g["series"] else "") + g["symbol"] for g in card["genes"][:GENES_SHOWN]]
        extra = len(card["genes"]) - GENES_SHOWN
        if genes:
            lines.append((tr("anatomy.viewer.genes", list=", ".join(genes) + (f" (+{extra})" if extra > 0 else "")), "normal"))
        if card["series_genes"]:
            lines.append((tr("anatomy.viewer.series", list=", ".join(card["series_genes"])), "gold"))
    lines.append((tr("anatomy.viewer.console"), "dim"))
    return {"title": shown_name(card), "subtitle": card["name"] if card.get("label") else "", "lines": lines}
