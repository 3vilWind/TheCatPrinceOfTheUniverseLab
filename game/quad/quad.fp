varying mediump vec2 var_texcoord0;

uniform lowp sampler2D QUAD_TEXTURE;
uniform lowp vec4 tint;
uniform lowp vec4 invert;

void main()
{
    // Pre-multiply alpha since all runtime textures already are
    //lowp vec4 tint_pm = vec4(tint.xyz * tint.w, tint.w);
    lowp vec4 textureColor = texture2D(QUAD_TEXTURE, var_texcoord0.xy);
    if (invert.r != 0.0) {
        gl_FragColor = invert - textureColor; // * tint_pm;
    } else {
        gl_FragColor = textureColor; // * tint_pm;
    }
}
