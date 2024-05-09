#version 330 core
uniform sampler2D inputTexture;
uniform int kernelSize;
uniform float sigma;
out vec4 FragColor;

float normpdf(in float x, in float sigma) {
    return 0.39894*exp(-0.5*x*x/(sigma*sigma))/sigma;
}

void main() {
    vec3 result = vec3(0.0);
    float w = 0.0;
    for(int i = -kernelSize; i <= kernelSize; ++i) {
        for(int j = -kernelSize; j <= kernelSize; ++j) {
            float weight = normpdf(length(vec2(float(i), float(j))), sigma);
            result += weight * texture(inputTexture, (gl_FragCoord.xy + vec2(float(i), float(j))) / textureSize(inputTexture, 0)).rgb;
            w += weight;
        }
    }
    FragColor = vec4(result / w, 1.0);
}
