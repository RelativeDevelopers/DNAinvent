"""Los textos del visor, dibujados con Pillow en una textura encima del cuerpo: la lista de capas
(a la izquierda; se pulsan para encender o apagar), lo que se ha elegido (a la derecha), la ayuda
y el aviso (abajo)."""

import numpy as np
from PIL import Image, ImageDraw

from microscope.hud import _font
from .scene import COLORS, hex_rgb

VERT = """#version 330
in vec2 in_pos; out vec2 v_uv;
void main(){ v_uv = vec2(in_pos.x*0.5+0.5, 0.5-in_pos.y*0.5); gl_Position = vec4(in_pos,0,1); }"""
FRAG = """#version 330
uniform sampler2D tex; in vec2 v_uv; out vec4 f;
void main(){ f = texture(tex, v_uv); }"""

WHITE, DIM, GOLD, PANEL = (235, 238, 242, 255), (150, 156, 166, 255), (255, 210, 70, 255), (12, 14, 18, 190)


def wrap(d, text: str, font, width: int) -> list[str]:
    """Parte el texto en líneas que caben en width píxeles con esa letra."""
    lines, current = [], ""
    for word in text.split():
        candidate = f"{current} {word}".strip()
        if current and d.textlength(candidate, font=font) > width:
            lines.append(current)
            current = word
        else:
            current = candidate
    return lines + ([current] if current else [])


class HUD:
    def __init__(self, ctx, tr):
        import moderngl
        self.ctx, self.tr = ctx, tr
        self.prog = ctx.program(vertex_shader=VERT, fragment_shader=FRAG)
        quad = np.array([-1, -1, 1, -1, -1, 1, 1, 1], dtype="f4")
        self.vao = ctx.vertex_array(self.prog, [(ctx.buffer(quad.tobytes()), "2f", "in_pos")])
        self._moderngl = moderngl
        self.tex = None
        self.show_help = True
        self.layer_rows: list[tuple[str, tuple[int, int, int, int]]] = []   # (sistema, caja en px)
        self._fonts = {}

    def fonts(self, k):
        key = round(k, 2)
        if key not in self._fonts:
            self._fonts[key] = (_font(int(24 * k)), _font(int(16 * k)), _font(int(13 * k)))
        return self._fonts[key]

    def layer_at(self, x: int, y_top: int) -> str | None:
        """La capa bajo el ratón (y contada desde arriba, como en la imagen)."""
        for system, (x0, y0, x1, y1) in self.layer_rows:
            if x0 <= x <= x1 and y0 <= y_top <= y1:
                return system
        return None

    def compose(self, w: int, h: int, layers: list[dict], info: dict | None, message: str = "") -> Image.Image:
        tr = self.tr
        k = max(h / 900, 0.75)
        big, font, small = self.fonts(k)
        img = Image.new("RGBA", (w, h), (0, 0, 0, 0))
        d = ImageDraw.Draw(img)
        pad, lh = int(14 * k), int(23 * k)

        # ---- capas
        self.layer_rows = []
        x, y = pad, pad
        width = int(250 * k)
        d.rounded_rectangle((x - 6, y - 6, x + width, y + lh * (len(layers) + 1) + 6), int(8 * k), fill=PANEL)
        d.text((x, y), tr("anatomy.viewer.layers"), font=small, fill=DIM)
        y += lh
        for layer in layers:
            on = layer["on"]
            color = tuple(int(c * 255) for c in hex_rgb(COLORS.get(layer["system"], "#9a9a9a"))) + (255 if on else 90,)
            box = int(13 * k)
            d.rectangle((x, y + 4, x + box, y + 4 + box), fill=color if on else None, outline=color, width=2)
            d.text((x + box + 10, y), tr(f"anatomy.systems.{layer['system']}"), font=font, fill=WHITE if on else DIM)
            self.layer_rows.append((layer["system"], (x - 6, y, x + width, y + lh)))
            y += lh

        # ---- lo elegido
        right_w = int(min(560 * k, w * 0.42))
        rx = w - right_w - pad
        lines = []
        if info is None:
            lines.append((font, DIM, tr("anatomy.viewer.none")))
        else:
            lines += [(big, WHITE, part) for part in wrap(d, info["title"], big, right_w - pad)]
            if info.get("subtitle"):
                lines.append((small, DIM, info["subtitle"]))
            for text, color in info["lines"]:
                for part in wrap(d, text, font, right_w - pad):
                    lines.append((font, GOLD if color == "gold" else DIM if color == "dim" else WHITE, part))
        height = sum(int(f.size * 1.35) for f, _, _ in lines) + 2 * pad
        d.rounded_rectangle((rx - pad, pad - 6, w - pad + 6, pad - 6 + height), int(8 * k), fill=PANEL)
        y = pad
        for f, color, text in lines:
            d.text((rx, y), text, font=f, fill=color)
            y += int(f.size * 1.35)

        # ---- ayuda, aviso y mensajes
        bottom = h - pad
        d.text((pad, bottom - small.size), tr("anatomy.viewer.footer"), font=small, fill=DIM)
        bottom -= int(small.size * 1.6)
        if self.show_help:
            for line in reversed(tr.lines("anatomy.viewer.help")):
                d.text((pad, bottom - font.size), line, font=font, fill=WHITE)
                bottom -= int(font.size * 1.4)
        if message:
            tw = d.textlength(message, font=big)
            d.rounded_rectangle(((w - tw) / 2 - pad, h * 0.08 - 6, (w + tw) / 2 + pad, h * 0.08 + big.size * 1.4),
                                int(8 * k), fill=PANEL)
            d.text(((w - tw) / 2, h * 0.08), message, font=big, fill=GOLD)
        return img

    def upload(self, img: Image.Image):
        if self.tex is None or self.tex.size != img.size:
            if self.tex is not None:
                self.tex.release()
            self.tex = self.ctx.texture(img.size, 4)
        self.tex.write(img.tobytes())

    def draw(self):
        if self.tex is None:
            return
        mgl = self._moderngl
        self.ctx.disable(mgl.DEPTH_TEST)
        self.ctx.enable(mgl.BLEND)
        self.ctx.blend_func = mgl.SRC_ALPHA, mgl.ONE_MINUS_SRC_ALPHA
        self.tex.use(0)
        self.vao.render(mgl.TRIANGLE_STRIP)
        self.ctx.disable(mgl.BLEND)
