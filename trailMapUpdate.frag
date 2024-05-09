#version 330 core
uniform sampler2D trailMap;
uniform float decay;
out vec4 FragColor;

void main() {
    vec2 uv = gl_FragCoord.xy / textureSize(trailMap, 0);
    float current = texture(trailMap, uv).r;
    FragColor = vec4(max(current - decay, 0.0), 0.0, 0.0, 1.0);
}
