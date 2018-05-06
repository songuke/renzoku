#ifndef _GL_VERTEX_H_
#define _GL_VERTEX_H_

#include "common.h"

namespace Renzoku {

struct GlVertex {
	Float pos[3];
	Float normal[3];
	Float kd[3];
	Float ks[3];	    
    Float tangent[3];
	Float glossy[2];  	    // exponent parameter for Ward BRDF		    
    Float material_type;	// Lambertian or Ward BRDF
    Float uv[2];            // all is float for easy loading
    Float padding[4];
};

} // end namespace

#endif
