#version 330
// Modo plástico: cada partícula es un cuadrado (dos triángulos) que mira siempre a la
// cámara, con la imagen de la bola de plástico. Se dibuja con instancias: un mismo
// cuadrado se repite en el centro de cada partícula.
in vec2 in_corner;      // esquina del cuadrado: (-1, -1) ... (1, 1)
in vec3 in_center;      // por instancia: centro de la bola (unidades del mundo)
in float in_radius;     // radio de la bola (unidades del mundo)
in float in_layer;      // qué imagen (capa de la textura): protón, neutrón, electrón...
in float in_alpha;
in float in_min_px;     // radio mínimo en píxeles (partículas puntuales: electrón, quarks)

uniform mat4 u_view;
uniform mat4 u_proj;
uniform float u_viewport_h;
uniform float u_ball;       // fracción del cuadrado que ocupa la bola en la imagen
uniform float u_near_fade;  // >0: la bola se vuelve transparente al acercarse (en radios)

out vec2 v_uv;
out float v_layer;
out float v_alpha;
out vec3 v_eye_center;
out float v_radius;

void main() {
    vec4 eye = u_view * vec4(in_center, 1.0);
    float dist = max(-eye.z, 1e-9);
    float px = in_radius * u_proj[1][1] * u_viewport_h * 0.5 / dist;
    // demasiado pequeña para verse: se agranda hasta in_min_px (como los puntos)
    float r = in_radius * max(1.0, in_min_px / max(px, 1e-6));
    float half_size = r / u_ball;
    gl_Position = u_proj * vec4(eye.xyz + vec3(in_corner * half_size, 0.0), 1.0);

    float a = in_alpha;
    if (u_near_fade > 0.0) {
        a *= clamp((dist / r - 1.2) / (u_near_fade - 1.2), 0.0, 1.0);
    }
    v_uv = in_corner;
    v_layer = in_layer;
    v_alpha = a;
    v_eye_center = eye.xyz;
    v_radius = r;
}
