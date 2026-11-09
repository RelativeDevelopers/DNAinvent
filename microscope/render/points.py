"""
Renderizado en GPU de nubes de puntos y de líneas (enlaces).

Todo se dibuja con mezcla aditiva: cada punto "emite luz" y donde se
acumulan muchos (zonas de alta densidad electrónica) se ve más brillante.
"""
from dataclasses import dataclass
from pathlib import Path

import moderngl
import numpy as np

from microscope.mathutil import gl

SHADERS = Path(__file__).resolve().parent / "shaders"


@dataclass
class Camera:
    view: np.ndarray
    proj: np.ndarray
    viewport_h: int
    time: float


@dataclass
class PointStyle:
    vib_amp: float = 0.0
    twinkle: float = 0.0
    min_px: float = 1.5
    max_px: float = 600.0
    softness: float = 1.0
    dim_small: bool = True     # False para partículas puntuales (quarks)


class Renderer:
    def __init__(self, ctx: moderngl.Context):
        self.ctx = ctx
        self.point_prog = ctx.program(
            vertex_shader=(SHADERS / "points.vert").read_text(),
            fragment_shader=(SHADERS / "points.frag").read_text())
        self.line_prog = ctx.program(
            vertex_shader=(SHADERS / "lines.vert").read_text(),
            fragment_shader=(SHADERS / "lines.frag").read_text())
        self._sprites = None          # modo plástico (render/sprites.py): se crea al usarlo

    # ------------------------------------------------------------ modo plástico
    def _sprite_renderer(self):
        if self._sprites is None:
            from microscope.render.sprites import SpriteRenderer
            self._sprites = SpriteRenderer(self.ctx)
        return self._sprites

    def sprite_program(self):
        return self._sprite_renderer().prog

    @property
    def sprite_quad(self):
        return self._sprite_renderer().quad

    def sprites(self, sprites, cam: "Camera", fade: float, near_fade: float = 0.0):
        """Bolas de plástico (SpriteSet). near_fade > 0: se vuelven transparentes al
        acercarse a menos de near_fade radios (para poder entrar en un protón)."""
        if sprites.count == 0 or fade <= 0.001:
            return
        self._sprite_renderer().draw(sprites, cam, fade, near_fade)

    def begin(self, clear=(0.004, 0.006, 0.016)):
        self.ctx.clear(*clear, 1.0)
        self.ctx.disable(moderngl.DEPTH_TEST)
        self.ctx.enable(moderngl.BLEND | moderngl.PROGRAM_POINT_SIZE)
        self.ctx.blend_func = moderngl.ONE, moderngl.ONE  # aditivo

    def points(self, cloud: "PointCloud", cam: Camera, style: PointStyle, fade: float,
               first: int = 0, count: int = -1):
        if cloud.count == 0 or fade <= 0.001:
            return
        p = self.point_prog
        p["u_view"].write(gl(cam.view))
        p["u_proj"].write(gl(cam.proj))
        p["u_viewport_h"].value = float(cam.viewport_h)
        p["u_time"].value = cam.time
        p["u_vib_amp"].value = style.vib_amp
        p["u_twinkle"].value = style.twinkle
        p["u_min_px"].value = style.min_px
        p["u_max_px"].value = style.max_px
        p["u_softness"].value = style.softness
        p["u_dim"].value = 1.0 if style.dim_small else 0.0
        p["u_fade"].value = fade
        cloud.vao.render(moderngl.POINTS, vertices=count, first=first)

    def lines(self, lines: "LineSet", cam: Camera, fade: float, dash_period=0.0):
        if lines.count == 0 or fade <= 0.001:
            return
        p = self.line_prog
        p["u_view"].write(gl(cam.view))
        p["u_proj"].write(gl(cam.proj))
        p["u_fade"].value = fade
        p["u_dash_period"].value = dash_period
        lines.vao.render(moderngl.LINES)


class PointCloud:
    """Puntos en GPU. Se crea una vez y se dibuja cada frame.
    Con dynamic=True se pueden reescribir cada frame con write()."""

    DTYPE = [("p", "f4", 3), ("c", "f4", 4), ("s", "f4"), ("k", "f4")]

    def __init__(self, renderer: Renderer, pos, color, size, seed, dynamic=False):
        n = len(pos)
        self.count = n
        if n == 0:
            return
        self.vbo = renderer.ctx.buffer(self._pack(pos, color, size, seed), dynamic=dynamic)
        self.vao = renderer.ctx.vertex_array(
            renderer.point_prog,
            [(self.vbo, "3f 4f 1f 1f", "in_pos", "in_color", "in_size", "in_seed")])

    def _pack(self, pos, color, size, seed):
        data = np.empty(len(pos), dtype=self.DTYPE)
        data["p"], data["c"], data["s"], data["k"] = pos, color, size, seed
        return data.tobytes()

    def write(self, pos, color, size, seed):
        """Actualiza todos los puntos (mismo número que al crearla)."""
        if self.count:
            self.vbo.write(self._pack(pos, color, size, seed))

    def release(self):
        if self.count:
            self.vao.release()
            self.vbo.release()


class LineSet:
    def __init__(self, renderer: Renderer, a, b, color_a, color_b=None):
        """a, b: extremos (N,3); colores (N,4)."""
        n = len(a)
        self.count = n
        if n == 0:
            return
        color_b = color_a if color_b is None else color_b
        length = np.linalg.norm(np.asarray(b) - np.asarray(a), axis=1)
        data = np.empty(2 * n, dtype=[("p", "f4", 3), ("c", "f4", 4), ("d", "f4")])
        data["p"][0::2], data["p"][1::2] = a, b
        data["c"][0::2], data["c"][1::2] = color_a, color_b
        data["d"][0::2], data["d"][1::2] = 0.0, length
        self.vbo = renderer.ctx.buffer(data.tobytes())
        self.vao = renderer.ctx.vertex_array(
            renderer.line_prog, [(self.vbo, "3f 4f 1f", "in_pos", "in_color", "in_dash")])

    def release(self):
        if self.count:
            self.vao.release()
            self.vbo.release()
