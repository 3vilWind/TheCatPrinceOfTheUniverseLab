varying mediump vec2 var_texcoord0;

uniform lowp sampler2D texture_sampler;
uniform lowp vec4 tint;
uniform lowp vec4 invert;

void main()
{
    // Pre-multiply alpha since all runtime textures already are
    lowp vec4 tint_pm = vec4(tint.xyz * tint.w, tint.w);
    lowp vec4 textureColor = texture2D(texture_sampler, var_texcoord0.xy);
    if (invert.x) {
        gl_FragColor = vec4(1.0 - textureColor.r, 1.0 - textureColor.g, 1.0 - textureColor.b, 1) * tint_pm;
    } else {
        gl_FragColor = textureColor * tint_pm;
    }
}
