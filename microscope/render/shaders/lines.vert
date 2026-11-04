#version 330
in vec3 in_pos;
in vec4 in_color;
in float in_dash;    // distancia recorrida a lo largo de la línea (para guiones)

uniform mat4 u_view;
uniform mat4 u_proj;

out vec4 v_color;
out float v_dash;

void main() {
    gl_Position = u_proj * u_view * vec4(in_pos, 1.0);
    v_color = in_color;
    v_dash = in_dash;
}
