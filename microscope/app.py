"""
La ventana del microscopio: controles, bucle y teclas.

La abre microscope/__main__.py (python -m microscope ATGC...) o la opción 7 de la consola.
"""
import moderngl
import numpy as np
import pyglet
from pyglet.window import key, mouse

from microscope.engine import Engine
from microscope.hud import HUD
from ui import settings


class App(pyglet.window.Window):
    def __init__(self, molecule, tr, didactic=True, comparison=None, label="", highlight=None,
                 mode=None):
        self.tr = tr
        # depth_size: el modo plástico usa el búfer de profundidad (las bolas se tapan)
        config = pyglet.gl.Config(double_buffer=True, major_version=3, minor_version=3,
                                  forward_compatible=True, sample_buffers=1, samples=4,
                                  depth_size=24)
        try:
            super().__init__(1600, 900, tr("microscope.title"), resizable=True, config=config)
        except pyglet.window.NoSuchConfigException:
            super().__init__(1600, 900, tr("microscope.title"), resizable=True)
        self.ctx = moderngl.create_context()
        print("GPU:", self.ctx.info["GL_RENDERER"])
        if "Intel" in self.ctx.info["GL_RENDERER"]:
            print(tr("microscope.console.integrated_gpu"))

        if getattr(molecule, "kind", "dna") == "protein":
            for mol in [molecule] + ([comparison.second] if comparison else []):
                print(tr("microscope.console.built_protein", residues=len(mol.sequence),
                         atoms=mol.n_atoms))
        else:
            print(tr("microscope.console.built", pairs=len(molecule.sequence),
                     atoms=molecule.n_atoms))
        self.engine = Engine(self.ctx, molecule, tr, didactic, comparison, label, highlight, mode)
        self.hud = HUD(self.ctx, tr)
        self.hud_timer = 0.0
        self.keys = key.KeyStateHandler()
        self.push_handlers(self.keys)
        self.captured = False
        self.hud.flash(tr("microscope.hud.click"), 0.0, 4.0)
        pyglet.clock.schedule_interval(self.tick, 1 / 144)

    # ------------------------------------------------------------ bucle
    def tick(self, dt):
        dt = min(dt, 0.05)
        k, tr, eng = self.keys, self.tr, self.engine
        move = np.array([
            (k[key.D] - k[key.A]),
            (k[key.SPACE] - (k[key.LCTRL] or k[key.RCTRL])),
            (k[key.W] - k[key.S]),
        ], dtype=float)
        manual = move.any() or k[key.Q] or k[key.E]
        ap = eng.autopilot
        if manual and ap.active:
            ap.cancel()
            self.hud.flash(tr("microscope.autopilot.off"), eng.time, 1.5)
        if k[key.Q]:
            eng.ship.roll(1.2 * dt)
        if k[key.E]:
            eng.ship.roll(-1.2 * dt)
        prev_level = eng.level
        eng.update(dt, move, k[key.LSHIFT] or k[key.RSHIFT])
        if eng.level is not prev_level:
            entering = eng.level.depth > prev_level.depth
            text = "microscope.hud.entering" if entering else "microscope.hud.leaving"
            self.hud.flash(tr(text, name=eng.level.name), eng.time)
        if ap.message:
            self.hud.flash(ap.message, eng.time, 3.0)
            ap.message = ""
        self.hud_timer -= dt
        if self.hud_timer <= 0:
            self.refresh_hud()

    def refresh_hud(self):
        self.hud_timer = 0.12
        w, h = self.get_framebuffer_size()
        self.hud.upload(self.hud.compose(self.engine.hud_state(w, h), w, h, self.engine.time))

    def on_draw(self):
        w, h = self.get_framebuffer_size()
        self.ctx.screen.use()
        self.engine.render(w, h)
        self.hud.draw()

    # ------------------------------------------------------------ entrada
    def on_mouse_press(self, x, y, button, modifiers):
        if button == mouse.LEFT and not self.captured:
            self.captured = True
            self.set_exclusive_mouse(True)

    def _manual_look(self, dx, dy):
        ap = self.engine.autopilot
        if ap.active and (abs(dx) + abs(dy)) > 2:      # ignorar temblores mínimos del ratón
            ap.cancel()
            self.hud.flash(self.tr("microscope.autopilot.off"), self.engine.time, 1.5)
        if not ap.active:
            self.engine.ship.look(dx, dy)

    def on_mouse_motion(self, x, y, dx, dy):
        if self.captured:
            self._manual_look(dx, dy)

    def on_mouse_drag(self, x, y, dx, dy, buttons, modifiers):
        self._manual_look(dx, dy)

    def on_mouse_scroll(self, x, y, scroll_x, scroll_y):
        s = self.engine.ship
        s.speed_mult = float(np.clip(s.speed_mult * (1.25 ** scroll_y), 0.05, 20))

    def on_key_press(self, symbol, modifiers):
        tr, eng = self.tr, self.engine
        if symbol == key.ESCAPE:
            if self.captured:
                self.captured = False
                self.set_exclusive_mouse(False)
            else:
                self.close()
            return pyglet.event.EVENT_HANDLED   # que Esc no cierre la ventana directamente
        if symbol in (key.ENTER, key.NUM_ENTER):
            eng.autopilot.dive(eng, full=bool(modifiers & key.MOD_SHIFT))
            if eng.autopilot.active:
                self.hud.flash(eng.autopilot.flash_text(), eng.time, 1.5)
        elif symbol == key.BACKSPACE:
            eng.autopilot.ascend(eng)
        elif symbol == key.TAB:
            mode = eng.toggle_mode(-1 if modifiers & key.MOD_SHIFT else 1)
            self.hud.flash(tr(f"microscope.hud.now_{mode}"), eng.time)
        elif symbol in (key.O, key.I):
            text = eng.cycle_subshell() if symbol == key.O else eng.cycle_orbital()
            if text is None:
                text = tr("microscope.hud.only_atomic" if symbol == key.O or eng.mode != "atomic"
                          else "microscope.hud.only_p")
            self.hud.flash(text, eng.time)
        elif symbol == key.C and eng.comparison is not None:
            if eng.toggle_layout():
                self.hud.flash(tr("microscope.hud.layout_overlay" if eng.layout == "overlay"
                                  else "microscope.hud.layout_split"), eng.time)
            else:
                self.hud.flash(tr("microscope.hud.only_at_protein"), eng.time)
        elif symbol == key.V and eng.comparison is not None:
            if eng.switch_side():
                self.hud.flash(tr("microscope.hud.side", side=eng.view.label), eng.time)
            else:
                self.hud.flash(tr("microscope.hud.only_at_protein"), eng.time)
        elif symbol == key.L:
            settings.save({"language": tr.next_language()})
            self.set_caption(tr("microscope.title"))
            self.hud.flash(tr("microscope.hud.language"), eng.time)
        elif symbol == key.H:
            self.hud.show_help = not self.hud.show_help
        elif symbol == key.F11:
            self.set_fullscreen(not self.fullscreen)
        self.refresh_hud()
