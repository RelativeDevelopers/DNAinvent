#version 330
// La imagen de la bola, con la profundidad de una esfera de verdad: así dos bolas que se
// tocan (los nucleones del núcleo) se cortan bien y la de delante tapa a la de detrás.
uniform sampler2DArray u_tex;
uniform mat4 u_proj;
uniform float u_ball;
uniform float u_fade;      // fundido entre niveles

in vec2 v_uv;
in float v_layer;
in float v_alpha;
in vec3 v_eye_center;
in float v_radius;
out vec4 f_color;

void main() {
    vec4 c = texture(u_tex, vec3(v_uv.x * 0.5 + 0.5, 0.5 - v_uv.y * 0.5, v_layer));
    float a = c.a * v_alpha * u_fade;
    if (a < 0.03) discard;

    vec2 q = v_uv / u_ball;                          // -1..1 sobre la bola
    float nz = sqrt(max(1.0 - dot(q, q), 0.0));
    vec3 p = v_eye_center + vec3(q, nz) * v_radius;  // punto de la superficie de la esfera
    vec4 clip = u_proj * vec4(p, 1.0);
    gl_FragDepth = clamp(clip.z / clip.w * 0.5 + 0.5, 0.0, 1.0);
    f_color = vec4(c.rgb, a);
}
