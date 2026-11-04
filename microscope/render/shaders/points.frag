#version 330
in vec4 v_color;
uniform float u_fade;      // fundido entre niveles
uniform float u_softness;  // 0 = disco nítido, 1 = brillo gaussiano
out vec4 f_color;

void main() {
    vec2 d = gl_PointCoord * 2.0 - 1.0;
    float r2 = dot(d, d);
    if (r2 > 1.0) discard;
    float glow = exp(-r2 * 4.0);
    float disk = smoothstep(1.0, 0.7, r2);
    float a = mix(disk, glow, u_softness) * v_color.a * u_fade;
    f_color = vec4(v_color.rgb * a, a);   // mezcla aditiva
}
