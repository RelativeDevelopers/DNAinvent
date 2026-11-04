#version 330
in vec4 v_color;
in float v_dash;
uniform float u_fade;
uniform float u_dash_period;   // 0 = línea continua
out vec4 f_color;

void main() {
    if (u_dash_period > 0.0 && fract(v_dash / u_dash_period) > 0.5) discard;
    float a = v_color.a * u_fade;
    f_color = vec4(v_color.rgb * a, a);
}
