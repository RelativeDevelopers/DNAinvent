"""El dibujo del cuerpo con moderngl, y saber qué pieza hay en un píxel.

Cada vértice lleva su número de pieza. Una textura de una fila (una columna por pieza) dice el
color de cada pieza y su estado: oculta, normal o elegida. Así encender una capa o elegir una
pieza solo cambia esa textura, no las mallas.

Para saber qué hay bajo el ratón se dibuja otra vez, fuera de la pantalla, con el número de
cada pieza como color, y se lee ese píxel.
"""

import numpy as np

HIDDEN, VISIBLE, SELECTED = 0, 128, 255
PLAIN = (0.86, 0.80, 0.74)               # «un solo color»: como una figura de escayola
HIGHLIGHT = (1.0, 0.82, 0.25)
BACKGROUND = (0.07, 0.08, 0.10)

VERT = """#version 330
uniform mat4 mvp; uniform mat4 view;
uniform sampler2D state;
in vec3 in_pos; in vec3 in_norm; in float in_piece;
out vec3 v_norm; out vec4 v_state;
void main(){
    v_state = texelFetch(state, ivec2(int(in_piece + 0.5), 0), 0);
    v_norm = mat3(view) * in_norm;
    gl_Position = (v_state.a < 0.25) ? vec4(2.0, 2.0, 2.0, 1.0) : mvp * vec4(in_pos, 1.0);
}"""
FRAG = """#version 330
uniform vec3 highlight; uniform float strength;
in vec3 v_norm; in vec4 v_state; out vec4 f;
void main(){
    if (v_state.a < 0.25) discard;
    vec3 n = normalize(v_norm);
    if (!gl_FrontFacing) n = -n;
    float key = max(dot(n, normalize(vec3(0.25, 0.35, 1.0))), 0.0);
    float rim = pow(1.0 - max(n.z, 0.0), 2.0);
    vec3 base = v_state.rgb;
    bool chosen = v_state.a > 0.9 && strength > 0.5;
    if (chosen) base = mix(base, highlight, 0.6);
    vec3 c = base * (0.28 + 0.72 * key) + vec3(0.10) * rim;
    if (chosen) c += highlight * 0.12;
    f = vec4(c, 1.0);
}"""
PICK_VERT = """#version 330
uniform mat4 mvp; uniform sampler2D state;
in vec3 in_pos; in float in_piece;
flat out int v_piece; out float v_alpha;
void main(){
    v_piece = int(in_piece + 0.5);
    v_alpha = texelFetch(state, ivec2(v_piece, 0), 0).a;
    gl_Position = (v_alpha < 0.25) ? vec4(2.0, 2.0, 2.0, 1.0) : mvp * vec4(in_pos, 1.0);
}"""
PICK_FRAG = """#version 330
flat in int v_piece; in float v_alpha; out vec4 f;
void main(){
    if (v_alpha < 0.25) discard;
    int n = v_piece + 1;
    f = vec4(float(n & 255), float((n >> 8) & 255), float((n >> 16) & 255), 255.0) / 255.0;
}"""


def encode_pick(rgb: bytes) -> int | None:
    n = rgb[0] | (rgb[1] << 8) | (rgb[2] << 16)
    return n - 1 if n else None


class Renderer:
    def __init__(self, ctx, scene):
        import moderngl
        self.ctx, self.scene = ctx, scene
        self.prog = ctx.program(vertex_shader=VERT, fragment_shader=FRAG)
        self.pick_prog = ctx.program(vertex_shader=PICK_VERT, fragment_shader=PICK_FRAG)
        self.vaos, self.pick_vaos = {}, {}
        for mesh in scene.meshes:
            vbo = ctx.buffer(np.hstack([mesh.positions, mesh.normals, mesh.piece[:, None]]).astype("f4").tobytes())
            ibo = ctx.buffer(mesh.indices.astype("u4").tobytes())
            self.vaos[mesh.system] = ctx.vertex_array(
                self.prog, [(vbo, "3f 3f 1f", "in_pos", "in_norm", "in_piece")], ibo)
            self.pick_vaos[mesh.system] = ctx.vertex_array(
                self.pick_prog, [(vbo, "3f 12x 1f", "in_pos", "in_piece")], ibo)
        n = len(scene.pieces)
        self.colors = np.array([_rgb(p["system"]) for p in scene.pieces], np.float32)
        self.state = np.zeros((n, 4), np.uint8)
        self.state_tex = ctx.texture((n, 1), 4, dtype="f1")
        self.state_tex.filter = (moderngl.NEAREST, moderngl.NEAREST)
        # la piel tapa todo lo demás: se empieza sin ella (tecla A o la lista de capas)
        self.visible_systems = set(scene.systems) - {"skin"}
        self.selected: set[int] = set()
        self.isolate = False
        self.plain = False
        self._pick_fbo, self._pick_size = None, None
        self.update_state()

    # ---------------------------------------------------------------- estado
    def update_state(self):
        colors = np.tile(np.array(PLAIN, np.float32), (len(self.colors), 1)) if self.plain else self.colors
        self.state[:, :3] = np.clip(colors * 255, 0, 255).astype(np.uint8)
        visible = np.array([p["system"] in self.visible_systems for p in self.scene.pieces])
        if self.isolate and self.selected:
            visible = np.zeros(len(visible), bool)
            visible[list(self.selected)] = True
        self.state[:, 3] = np.where(visible, VISIBLE, HIDDEN)
        if self.selected:
            chosen = np.zeros(len(visible), bool)
            chosen[list(self.selected)] = True
            self.state[chosen & visible, 3] = SELECTED
        self.state_tex.write(self.state.tobytes())
        systems = np.array([p["system"] for p in self.scene.pieces])
        self.drawn = set(systems[self.state[:, 3] != HIDDEN])

    def select(self, numbers):
        self.selected = set(numbers)
        if not self.selected:
            self.isolate = False
        self.update_state()

    def toggle_system(self, system: str):
        self.visible_systems ^= {system}
        self.update_state()

    def show_all(self):
        self.visible_systems = set(self.scene.systems)
        self.isolate = False
        self.update_state()

    # ---------------------------------------------------------------- dibujo
    def _matrices(self, camera, w, h):
        view = camera.view()
        proj = camera.projection(w / max(h, 1))
        return (proj @ view).T.astype("f4"), view.T.astype("f4")

    def render(self, camera, w, h, fbo=None):
        import moderngl
        target = fbo or self.ctx.screen
        target.use()
        self.ctx.viewport = (0, 0, w, h)
        self.ctx.clear(*BACKGROUND, 1.0)
        self.ctx.enable(moderngl.DEPTH_TEST)
        self.ctx.disable(moderngl.CULL_FACE)          # algunas mallas vienen con caras al revés
        mvp, view = self._matrices(camera, w, h)
        self.prog["mvp"].write(mvp.tobytes())
        self.prog["view"].write(view.tobytes())
        # al ver solo lo elegido no hace falta resaltarlo: se ven sus colores
        self.prog["highlight"].value = HIGHLIGHT
        self.prog["strength"].value = 0.0 if self.isolate else 1.0
        self.state_tex.use(0)
        self.prog["state"].value = 0
        for system, vao in self.vaos.items():
            if system in self.drawn:
                vao.render()

    def pick(self, camera, w, h, x, y) -> int | None:
        """Número de la pieza en el píxel (x, y), con y desde abajo, o None."""
        import moderngl
        if self._pick_size != (w, h):
            self._pick_fbo = self.ctx.simple_framebuffer((w, h), components=4)
            self._pick_size = (w, h)
        previous = self.ctx.fbo
        fbo = self._pick_fbo
        fbo.use()
        self.ctx.viewport = (0, 0, w, h)
        fbo.clear(0, 0, 0, 0, depth=1.0)
        self.ctx.enable(moderngl.DEPTH_TEST)
        self.ctx.disable(moderngl.CULL_FACE)
        mvp, _ = self._matrices(camera, w, h)
        self.pick_prog["mvp"].write(mvp.tobytes())
        self.state_tex.use(0)
        self.pick_prog["state"].value = 0
        for system, vao in self.pick_vaos.items():
            if system in self.drawn:
                vao.render()
        x, y = int(np.clip(x, 0, w - 1)), int(np.clip(y, 0, h - 1))
        number = encode_pick(fbo.read(viewport=(x, y, 1, 1), components=3))
        if previous is not None:
            previous.use()
        return number if number is not None and number < len(self.scene.pieces) else None


def _rgb(system: str):
    from .scene import COLORS, hex_rgb
    return hex_rgb(COLORS.get(system, "#9a9a9a"))
