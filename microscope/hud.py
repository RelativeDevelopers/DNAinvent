"""
HUD: información superpuesta, dibujada con Pillow en una textura.

Se redibuja unas 8 veces por segundo (no hace falta más para texto). Todos los textos
salen del traductor (i18n/locales/*.json), así que cambian al momento con la tecla L.
"""
import math

import moderngl
import numpy as np
from PIL import Image, ImageDraw, ImageFont

from microscope import scale

FONT_CANDIDATES = ["consola.ttf", "CascadiaMono.ttf", "DejaVuSansMono.ttf",
                   "LiberationMono-Regular.ttf", "Menlo.ttc"]

DEPTH_MARKS = [(2e-9, "dna"), (1e-10, "atom"), (3e-15, "nucleus"), (8.4e-16, "proton"),
               (1e-18, "quarks")]
FIRST_MARK = {"dna": (2e-9, "dna"), "protein": (5e-9, "protein")}   # según la molécula
BAR_TOP, BAR_BOTTOM = 1e-6, 1e-19
SEQUENCE_SHOWN = 40          # bases de la secuencia que caben arriba

VERT = """#version 330
in vec2 in_pos; out vec2 v_uv;
void main(){ v_uv = vec2(in_pos.x*0.5+0.5, 0.5-in_pos.y*0.5); gl_Position = vec4(in_pos,0,1); }"""
FRAG = """#version 330
uniform sampler2D tex; in vec2 v_uv; out vec4 f;
void main(){ f = texture(tex, v_uv); }"""


def _font(size):
    for name in FONT_CANDIDATES:
        try:
            return ImageFont.truetype(name, size)
        except OSError:
            continue
    return ImageFont.load_default(size=size)


class HUD:
    def __init__(self, ctx: moderngl.Context, tr):
        self.ctx = ctx
        self.tr = tr
        self.prog = ctx.program(vertex_shader=VERT, fragment_shader=FRAG)
        quad = np.array([-1, -1, 1, -1, -1, 1, 1, 1], dtype="f4")
        self.vbo = ctx.buffer(quad.tobytes())
        self.vao = ctx.vertex_array(self.prog, [(self.vbo, "2f", "in_pos")])
        self.tex = None
        self.size = (0, 0)
        self.show_help = True
        self.message = ""
        self.message_until = 0.0
        self._fonts = {}

    def flash(self, text, now, seconds=2.5):
        self.message, self.message_until = text, now + seconds

    def fonts(self, k):
        key = round(k, 2)
        if key not in self._fonts:
            self._fonts[key] = (_font(int(22 * k)), _font(int(16 * k)), _font(int(13 * k)))
        return self._fonts[key]

    # ------------------------------------------------------------ composición
    def compose(self, st: dict, w: int, h: int, now: float) -> Image.Image:
        tr = self.tr
        dec = tr.decimal_separator
        k = max(h / 900, 0.75)
        f_big, f, f_small = self.fonts(k)
        img = Image.new("RGBA", (w, h), (0, 0, 0, 0))
        d = ImageDraw.Draw(img)
        pad = int(16 * k)
        lh = int(21 * k)

        # --- panel superior izquierdo
        lvl = st["level"]
        lines = [
            (f_big, tr("microscope.hud.level", n=lvl.depth + 1, name=lvl.name), (255, 255, 255, 255)),
            (f, tr("microscope.hud.visible", width=scale.length(st["visible_width_m"], dec),
                   sci=tr.sci(st["visible_width_m"])), (190, 215, 255, 255)),
            (f, tr("microscope.hud.distance", d=scale.length(st["depth_m"], dec)), (190, 215, 255, 255)),
            (f, tr("microscope.hud.speed", v=scale.length(st["speed_m_s"], dec),
                   mult=f"{st['speed_mult']:.2g}".replace(".", dec)), (190, 215, 255, 255)),
            (f_small, tr("microscope.hud.time", t=st["time"]), (150, 170, 200, 255)),
            (f_small, tr("microscope.hud.mode_didactic" if st["didactic"] else "microscope.hud.mode_realistic"),
             (255, 210, 120, 255) if st["didactic"] else (140, 230, 200, 255)),
        ]
        y = pad
        box_w = max(int(560 * k), int(max(d.textlength(t, font=fo) for fo, t, _ in lines)) + pad)
        d.rounded_rectangle([pad // 2, pad // 2, pad + box_w, pad + lh * len(lines) + lh // 2],
                            radius=int(8 * k), fill=(5, 10, 25, 150))
        for font, text, col in lines:
            d.text((pad, y), text, font=font, fill=col)
            y += lh + (6 if font is f_big else 0)

        # --- objeto delante
        foc = st["focus"]
        if foc:
            y += lh // 2
            fl = [(f, "▶ " + foc.title, (255, 255, 255, 255))]
            if foc.detail:
                fl.append((f_small, foc.detail, (200, 210, 230, 255)))
            for e in foc.explanation:
                fl.append((f_small, e, (255, 220, 140, 255)))
            fw = max(int(620 * k), int(max(d.textlength(t, font=fo) for fo, t, _ in fl)) + pad)
            d.rounded_rectangle([pad // 2, y - pad // 3, pad + fw, y + lh * len(fl) + 4],
                                radius=int(8 * k), fill=(5, 10, 25, 150))
            for font, text, col in fl:
                d.text((pad, y), text, font=font, fill=col)
                y += lh

        # --- secuencia arriba en el centro, con el par (o aminoácido) actual resaltado
        mol = st["molecule"]
        self._sequence(d, mol, foc.pair_index if foc else None, w, pad, f, f_small)

        # --- pantalla partida: separación y nombre de cada mitad
        panels = st.get("panels") or []
        for panel in panels:
            self._panel(d, panel, w, h, pad, lh, f, k)

        # --- barra de profundidad a la derecha
        self._depth_bar(d, st["depth_m"], w, h, k, f_small, getattr(mol, "kind", "dna"))

        # --- mira central (una en cada mitad si la pantalla está partida)
        centers = ([(int((p["x0"] + p["x1"]) / 2 * w), h // 2) for p in panels] if panels
                   else [(w // 2, h // 2)])
        r = int(9 * k)
        for cx, cy in centers:
            for a, b in (((cx - r, cy), (cx - 3, cy)), ((cx + 3, cy), (cx + r, cy)),
                         ((cx, cy - r), (cx, cy - 3)), ((cx, cy + 3), (cx, cy + r))):
                d.line([a, b], fill=(255, 255, 255, 140))

        # --- piloto automático: estado y marcador sobre el objetivo
        if st.get("autopilot"):
            txt = "◆ " + st["autopilot"]
            tw = d.textlength(txt, font=f)
            d.rounded_rectangle([(w - tw) / 2 - 12, h - pad - lh * 4.6, (w + tw) / 2 + 12,
                                 h - pad - lh * 3.4], radius=int(6 * k), fill=(40, 30, 5, 170))
            d.text(((w - tw) / 2, h - pad - lh * 4.45), txt, font=f, fill=(255, 215, 110, 255))
        ts = st.get("target_screen")
        if ts and 0 <= ts[0] <= w and 0 <= ts[1] <= h:
            x, y, s_ = ts[0], ts[1], int(16 * k)
            col = (255, 215, 110, 230)
            for sx, sy in ((-1, -1), (1, -1), (-1, 1), (1, 1)):
                cx_, cy_ = x + sx * s_, y + sy * s_
                d.line([cx_, cy_, cx_ - sx * s_ // 2, cy_], fill=col, width=2)
                d.line([cx_, cy_, cx_, cy_ - sy * s_ // 2], fill=col, width=2)

        if self.show_help:
            help_lines = tr.lines("microscope.hud.help")
            if st.get("comparing"):
                help_lines = help_lines + [tr("microscope.hud.help_compare")]
            for i, text in enumerate(help_lines):
                d.text((pad, h - pad - lh * (len(help_lines) - i)), text, font=f_small,
                       fill=(160, 175, 200, 220))
        if now < self.message_until:
            tw = d.textlength(self.message, font=f_big)
            d.text(((w - tw) / 2, h * 0.72), self.message, font=f_big, fill=(255, 255, 255, 230))
        return img

    def _sequence(self, d, mol, pair, w, pad, font, f_small):
        """La secuencia (o el tramo), con la base del par que tienes delante en amarillo.

        pair es la posición en la secuencia completa; mol.start dice dónde empieza el tramo."""
        seq, first = mol.sequence, mol.start
        local = None if pair is None else pair - first
        start = 0
        if local is not None and len(seq) > SEQUENCE_SHOWN:
            start = int(np.clip(local - SEQUENCE_SHOWN // 2, 0, len(seq) - SEQUENCE_SHOWN))
        shown = seq[start:start + SEQUENCE_SHOWN]
        cw = d.textlength("M", font=font)
        x0 = (w - cw * len(shown)) / 2
        for i, ch in enumerate(shown):
            hi = local is not None and start + i == local
            d.text((x0 + i * cw, pad), ch, font=font,
                   fill=(255, 230, 120, 255) if hi else (120, 140, 170, 200))
            if hi:
                d.line([x0 + i * cw, pad + cw * 1.8, x0 + (i + 1) * cw, pad + cw * 1.8],
                       fill=(255, 230, 120, 255), width=2)
        # posición: "par 7 de 14" o, en un tramo, "par 127 · tramo 121–180 de 5.000"
        # (en una proteína: "aminoácido 7 de 147")
        tr = self.tr
        unit = "residue" if getattr(mol, "kind", "dna") == "protein" else "pair"
        if pair is not None:
            label = tr(f"microscope.hud.{unit}", pair=pair + 1, length=mol.length)
        else:
            label = ""
        if first > 0 or len(seq) < mol.length:
            span = tr("microscope.hud.span", first=first + 1, last=mol.end, length=mol.length)
            label = f"{label} · {span}" if label else span
        if label:
            tw = d.textlength(label, font=f_small)
            d.text(((w - tw) / 2, pad + cw * 2.3), label, font=f_small, fill=(150, 170, 200, 220))

    def _panel(self, d, panel, w, h, pad, lh, font, k):
        """Nombre de una mitad de la pantalla partida; la que manejas, resaltada."""
        x0, x1 = int(panel["x0"] * w), int(panel["x1"] * w)
        if panel["x0"] > 0:
            d.line([x0, 0, x0, h], fill=(120, 140, 170, 160), width=2)
        text = ("▶ " if panel["active"] else "") + panel["label"]
        tw = d.textlength(text, font=font)
        y = int(h * 0.80)
        col = (255, 230, 120, 255) if panel["active"] else (170, 185, 210, 220)
        x = (x0 + x1 - tw) / 2
        d.rounded_rectangle([x - 10, y - 4, x + tw + 10, y + lh], radius=int(6 * k),
                            fill=(5, 10, 25, 170))
        d.text((x, y), text, font=font, fill=col)

    def _depth_bar(self, d, depth_m, w, h, k, font, kind="dna"):
        x = w - int(40 * k)
        top, bot = int(80 * k), h - int(80 * k)
        lt, lb = math.log10(BAR_TOP), math.log10(BAR_BOTTOM)

        def y_of(m):
            return top + (lt - math.log10(m)) / (lt - lb) * (bot - top)

        d.line([x, top, x, bot], fill=(120, 140, 170, 160), width=2)
        for m, key in [FIRST_MARK[kind]] + DEPTH_MARKS[1:]:
            label = self.tr(f"microscope.depth_marks.{key}")
            y = y_of(m)
            d.line([x - 6, y, x + 6, y], fill=(150, 170, 200, 200), width=2)
            tw = d.textlength(label, font=font)
            d.text((x - 14 - tw, y - 8 * k), label, font=font, fill=(150, 170, 200, 220))
        y = y_of(min(max(depth_m, BAR_BOTTOM), BAR_TOP))
        d.polygon([(x + 8, y), (x + 20, y - 7), (x + 20, y + 7)], fill=(255, 220, 120, 255))

    # ------------------------------------------------------------ GPU
    def upload(self, img: Image.Image):
        if self.tex is None or self.size != img.size:
            if self.tex:
                self.tex.release()
            self.tex = self.ctx.texture(img.size, 4)
            self.size = img.size
        self.tex.write(img.tobytes())

    def draw(self):
        if not self.tex:
            return
        self.ctx.enable(moderngl.BLEND)
        self.ctx.blend_func = moderngl.SRC_ALPHA, moderngl.ONE_MINUS_SRC_ALPHA
        self.tex.use(0)
        self.vao.render(moderngl.TRIANGLE_STRIP)
