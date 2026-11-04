#version 330
// Cada punto: posición, color+brillo, tamaño en unidades del mundo y una
// "semilla" para animar (vibración térmica del átomo y parpadeo).
in vec3 in_pos;
in vec4 in_color;
in float in_size;
in float in_seed;

uniform mat4 u_view;
uniform mat4 u_proj;
uniform float u_viewport_h;
uniform float u_time;
uniform float u_vib_amp;      // amplitud de vibración (unidades del mundo)
uniform float u_twinkle;      // 0..1: parpadeo estadístico de las nubes
uniform float u_min_px;
uniform float u_max_px;
uniform float u_dim;          // 1: atenuar puntos subpíxel; 0: brillo fijo (partículas puntuales)

out vec4 v_color;

float hash(float n) { return fract(sin(n) * 43758.5453123); }

void main() {
    vec3 p = in_pos;
    if (u_vib_amp > 0.0) {
        // tres osciladores con frecuencias distintas por átomo
        vec3 ph = vec3(hash(in_seed), hash(in_seed + 1.3), hash(in_seed + 2.7)) * 6.2831;
        vec3 fr = vec3(5.0, 6.3, 7.1) * (0.8 + 0.4 * hash(in_seed + 4.1));
        p += u_vib_amp * sin(u_time * fr + ph);
    }
    vec4 eye = u_view * vec4(p, 1.0);
    gl_Position = u_proj * eye;

    float dist = max(-eye.z, 1e-9);
    float px = in_size * u_proj[1][1] * u_viewport_h * 0.5 / dist;
    gl_PointSize = clamp(px, u_min_px, u_max_px);

    // puntos más pequeños que el mínimo: se atenúan en vez de agrandarse,
    // así la cantidad de "luz" se conserva al alejarse
    float dim = mix(1.0, clamp(px / u_min_px, 0.0, 1.0), u_dim);
    float tw = 1.0;
    if (u_twinkle > 0.0) {
        float step_t = floor(u_time * 12.0);
        tw = mix(1.0, 0.25 + 1.5 * hash(in_seed * 17.0 + step_t), u_twinkle);
    }
    v_color = vec4(in_color.rgb, in_color.a * dim * tw);
}
