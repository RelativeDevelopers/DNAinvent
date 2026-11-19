"""La ventana del visor de «Cuerpo humano» (pyglet + moderngl)."""

import json
import time

import moderngl
import pyglet
from pyglet.window import key, mouse

from dnainvent import anatomy
from dnainvent.anatomy import cards
from ui import settings

from .camera import OrbitCamera
from .hud import HUD
from .info import selection_info
from .renderer import Renderer

CLICK_PIXELS = 4             # si el ratón se mueve menos que esto, es un clic (no un giro)


class App(pyglet.window.Window):
    def __init__(self, scene, conn, tr, series=frozenset(), select: str | None = None):
        self.tr, self.scene, self.conn, self.series = tr, scene, conn, series
        config = pyglet.gl.Config(double_buffer=True, major_version=3, minor_version=3,
                                  forward_compatible=True, sample_buffers=1, samples=4, depth_size=24)
        try:
            super().__init__(1600, 900, tr("anatomy.viewer.title"), resizable=True, config=config)
        except pyglet.window.NoSuchConfigException:
            super().__init__(1600, 900, tr("anatomy.viewer.title"), resizable=True,
                             config=pyglet.gl.Config(double_buffer=True, depth_size=24))
        self.ctx = moderngl.create_context()
        self.renderer = Renderer(self.ctx, scene)
        self.camera = OrbitCamera(scene.radius)
        self.hud = HUD(self.ctx, tr)
        self.concept = None                  # FMA de lo elegido
        self.info = None
        self.message, self.message_until = "", 0.0
        self.press = None
        self.dirty = True
        print(tr("anatomy.viewer.loaded", pieces=len(scene.pieces), triangles=scene.triangles))
        if select:                           # desde la consola: se ve sola; con I, en su sitio
            self.renderer.isolate = True
            self.choose(select, frame=True)
        pyglet.clock.schedule_interval(self.tick, 1 / 30)

    # ---------------------------------------------------------------- elegir
    def choose(self, fma: str | None, frame: bool = False, piece: str | None = None):
        self.concept = fma
        if fma is None:
            self.renderer.select([])
            self.info = None
        else:
            ids = cards.concept_pieces(self.conn, fma) or ([piece] if piece else [])
            numbers = [self.scene.by_id[i] for i in ids if i in self.scene.by_id]
            self.renderer.select(numbers)
            for n in numbers:                # que se vea aunque su capa estuviera apagada
                self.renderer.visible_systems.add(self.scene.pieces[n]["system"])
            self.renderer.update_state()
            if frame and numbers:
                self.camera.frame(*self.scene.bounds(numbers))
            self.refresh_info()
            save_selection(fma, piece)
        self.dirty = True

    def refresh_info(self):
        card = cards.concept_card(self.conn, self.concept, self.tr.language, self.series) if self.concept else None
        self.info = selection_info(card, self.tr)

    def flash(self, text: str, seconds: float = 2.0):
        self.message, self.message_until = text, time.time() + seconds
        self.dirty = True

    # ---------------------------------------------------------------- dibujo
    def tick(self, dt):
        if self.message and time.time() > self.message_until:
            self.message = ""
            self.dirty = True

    def on_draw(self):
        w, h = self.get_framebuffer_size()
        self.renderer.render(self.camera, w, h)
        if self.dirty:
            layers = [{"system": s, "on": s in self.renderer.visible_systems} for s in self.scene.systems]
            self.hud.upload(self.hud.compose(w, h, layers, self.info, self.message))
            self.dirty = False
        self.hud.draw()

    def on_resize(self, width, height):
        self.dirty = True
        return super().on_resize(width, height)

    # ---------------------------------------------------------------- ratón
    def _scale(self):
        """Píxeles de la ventana -> píxeles del framebuffer (pantallas de alta densidad)."""
        w, h = self.get_framebuffer_size()
        return w / max(self.width, 1), h / max(self.height, 1)

    def on_mouse_press(self, x, y, button, modifiers):
        self.press = (x, y, button)

    def on_mouse_drag(self, x, y, dx, dy, buttons, modifiers):
        if buttons & mouse.LEFT and not (modifiers & key.MOD_SHIFT):
            self.camera.orbit(dx, -dy)
        else:
            sx, sy = self._scale()
            self.camera.pan(dx * sx, dy * sy, self.get_framebuffer_size()[1])

    def on_mouse_release(self, x, y, button, modifiers):
        if self.press is None or button != mouse.LEFT:
            return
        px, py, _ = self.press
        self.press = None
        if abs(x - px) + abs(y - py) > CLICK_PIXELS:
            return
        sx, sy = self._scale()
        w, h = self.get_framebuffer_size()
        fx, fy = int(x * sx), int(y * sy)
        system = self.hud.layer_at(fx, h - fy)
        if system:
            self.renderer.toggle_system(system)
            self.dirty = True
            return
        number = self.renderer.pick(self.camera, w, h, fx, fy)
        if number is None:
            self.choose(None)
        else:
            piece = self.scene.pieces[number]
            self.choose(piece["fma"], piece=piece["id"])

    def on_mouse_scroll(self, x, y, scroll_x, scroll_y):
        self.camera.zoom(scroll_y)

    # ---------------------------------------------------------------- teclas
    def on_key_press(self, symbol, modifiers):
        tr, r = self.tr, self.renderer
        if symbol == key.ESCAPE:
            if self.concept:
                self.choose(None)
            else:
                self.close()
            return pyglet.event.EVENT_HANDLED
        if symbol == key.F:
            if r.selected:
                self.camera.frame(*self.scene.bounds(r.selected))
            else:
                self.flash(tr("anatomy.viewer.nothing_selected"))
        elif symbol == key.I:
            if r.selected:
                r.isolate = not r.isolate
                r.update_state()
                if r.isolate:
                    self.camera.frame(*self.scene.bounds(r.selected))
                self.flash(tr("anatomy.viewer.isolate_on" if r.isolate else "anatomy.viewer.isolate_off"))
            else:
                self.flash(tr("anatomy.viewer.nothing_selected"))
        elif symbol == key.C:
            r.plain = not r.plain
            r.update_state()
            self.flash(tr("anatomy.viewer.colors_plain" if r.plain else "anatomy.viewer.colors_systems"))
        elif symbol == key.A:
            r.show_all()
            self.flash(tr("anatomy.viewer.all_layers"))
        elif symbol == key.R:
            self.camera.reset()
            self.flash(tr("anatomy.viewer.reset"))
        elif symbol == key.L:
            settings.save({"language": tr.next_language()})
            self.set_caption(tr("anatomy.viewer.title"))
            self.refresh_info()
            self.flash(tr("anatomy.viewer.language"))
        elif symbol == key.H:
            self.hud.show_help = not self.hud.show_help
        elif symbol == key.F11:
            self.set_fullscreen(not self.fullscreen)
        self.dirty = True


def save_selection(fma: str, piece: str | None):
    """Para que la consola enseñe la ficha de lo último elegido (Cuerpo humano → 4)."""
    try:
        anatomy.selection_file().write_text(json.dumps({"fma": fma, "piece": piece, "time": time.time()}),
                                            encoding="utf-8")
    except OSError:
        pass
