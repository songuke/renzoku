uniform sampler2D tex_output;

void main()
{	
	vec3 color = texture2D(tex_output, gl_TexCoord[0].xy).rgb;
	    
	gl_FragColor.rgb = pow(max(color, vec3(0.0)), vec3(1.0 / 2.2)); // TODO: perform proper normalization (with min/max) here	
	gl_FragColor.w = 1.0;
}
