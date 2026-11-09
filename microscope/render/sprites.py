"""
Modo plástico: partículas dibujadas como bolas de plástico brillante con su símbolo.

    protón  (+, rojo)      neutrón (0, azul grisáceo)     electrón (e⁻, verde)
    quark up (u, azul)     quark down (d, azul)           gluón (g, naranja)

Las imágenes están en render/textures/ (256 × 256, fondo transparente). Cada partícula es
un cuadrado que mira a la cámara con su imagen encima ("billboard"), pero con la
profundidad de una esfera: las bolas que se solapan se cortan como bolas de verdad.

Se usan cuadrados con instancias y no puntos (GL_POINTS) porque muchas tarjetas limitan
el tamaño de un punto a 64 o 256 píxeles, y al acercarse a un protón la bola ocupa toda
la pantalla.

Es un modo de juguete: ni los electrones ni los quarks tienen tamaño (son puntuales) ni
color. Los textos del HUD lo dicen.
"""
from pathlib import Path

import moderngl
import numpy as np
from PIL import Image

from microscope.mathutil import gl

SHADERS = Path(__file__).resolve().parent / "shaders"
TEXTURES = Path(__file__).resolve().parent / "textures"
KINDS = ("proton", "neutron", "electron", "up", "down", "gluon")
LAYER = {name: float(i) for i, name in enumerate(KINDS)}
BALL_FRACTION = 0.92          # en las imágenes, la bola ocupa el 92 % del ancho
TEXTURE_SIZE = 256


def texture_paths() -> list[Path]:
    return [TEXTURES / f"{name}.png" for name in KINDS]


class SpriteSet:
    """Un grupo de bolas. Se reescribe entero cuando cambia (cada frame si se mueven).

    pos (N,3), radius (N,) en unidades del mundo, kind: nombres de KINDS, alpha (N,),
    min_px (N,): radio mínimo en pantalla (los electrones y los quarks no tienen tamaño:
    sin mínimo, desaparecerían)."""

    DTYPE = [("p", "f4", 3), ("r", "f4"), ("k", "f4"), ("a", "f4"), ("m", "f4")]

    def __init__(self, renderer, pos, radius, kind, alpha=1.0, min_px=0.0):
        self.renderer = renderer
        self.count = len(pos)
        self.data = np.zeros(self.count, dtype=self.DTYPE)
        self.data["k"] = [LAYER[k] for k in kind]
        self.data["m"] = min_px
        self.set(pos, radius, alpha)
        self.vbo = self.vao = None
        if self.count:
            prog = renderer.sprite_program()
            self.vbo = renderer.ctx.buffer(reserve=self.data.nbytes, dynamic=True)
            self.vao = renderer.ctx.vertex_array(prog, [
                (renderer.sprite_quad, "2f", "in_corner"),
                (self.vbo, "3f 1f 1f 1f 1f/i", "in_center", "in_radius", "in_layer", "in_alpha",
                 "in_min_px")])

    def set(self, pos=None, radius=None, alpha=None):
        if pos is not None:
            self.data["p"] = pos
        if radius is not None:
            self.data["r"] = radius
        if alpha is not None:
            self.data["a"] = alpha

    @property
    def positions(self) -> np.ndarray:
        return self.data["p"]

    def release(self):
        if self.vao is not None:
            self.vao.release()
            self.vbo.release()
            self.vao = self.vbo = None


class SpriteRenderer:
    """Lo que el Renderer necesita para las bolas. Se crea la primera vez que se usa: en
    los modos didáctico y realista no se carga ninguna imagen."""

    def __init__(self, ctx: moderngl.Context):
        self.ctx = ctx
        self.prog = ctx.program(vertex_shader=(SHADERS / "sprites.vert").read_text(),
                                fragment_shader=(SHADERS / "sprites.frag").read_text())
        corners = np.array([-1, -1, 1, -1, -1, 1, 1, 1], dtype="f4")
        self.quad = ctx.buffer(corners.tobytes())
        layers = [Image.open(p).convert("RGBA").resize((TEXTURE_SIZE, TEXTURE_SIZE), Image.LANCZOS)
                  for p in texture_paths()]
        data = b"".join(im.tobytes() for im in layers)
        self.texture = ctx.texture_array((TEXTURE_SIZE, TEXTURE_SIZE, len(layers)), 4, data)
        self.texture.build_mipmaps()
        self.texture.filter = moderngl.LINEAR_MIPMAP_LINEAR, moderngl.LINEAR
        self.texture.repeat_x = self.texture.repeat_y = False

    def draw(self, sprites: SpriteSet, cam, fade: float, near_fade: float):
        ctx, p = self.ctx, self.prog
        # de atrás hacia delante: así la transparencia de los bordes (y del fundido) sale bien
        eye_z = sprites.data["p"] @ cam.view[2, :3] + cam.view[2, 3]
        sprites.vbo.write(sprites.data[np.argsort(eye_z)].tobytes())
        p["u_view"].write(gl(cam.view))
        p["u_proj"].write(gl(cam.proj))
        p["u_viewport_h"].value = float(cam.viewport_h)
        p["u_ball"].value = BALL_FRACTION
        p["u_near_fade"].value = float(near_fade)
        p["u_fade"].value = float(fade)
        p["u_tex"].value = 0
        self.texture.use(0)

        _clear_depth(ctx)               # cada nivel tiene sus unidades: profundidad nueva
        ctx.enable(moderngl.DEPTH_TEST)
        ctx.blend_func = moderngl.SRC_ALPHA, moderngl.ONE_MINUS_SRC_ALPHA
        sprites.vao.render(moderngl.TRIANGLE_STRIP, instances=sprites.count)
        ctx.disable(moderngl.DEPTH_TEST)
        ctx.blend_func = moderngl.ONE, moderngl.ONE          # vuelve la mezcla aditiva


def _clear_depth(ctx):
    """Borra solo la profundidad (y solo en el viewport: con pantalla partida, cada mitad)."""
    fbo = ctx.fbo
    viewport = ctx.viewport
    mask = fbo.color_mask
    fbo.color_mask = (False, False, False, False)
    fbo.clear(depth=1.0, viewport=viewport)
    fbo.color_mask = mask
    fbo.use()                  # vuelve a aplicar la máscara de color (si no, no se dibuja nada)
    ctx.viewport = viewport    # use() devuelve el viewport a la ventana entera
