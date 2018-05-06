#ifndef _GL_VIEW_H_
#define _GL_VIEW_H_

#include "view.h"

namespace Renzoku {

/**
 * OpenGL specific view. It can be used for general setups such as
 * frame buffer objects (FBO).
 * 
 * Currently leave as blank.
 */
class GLView : public View {
public:
    GLView(Scene *scene) : View(scene) {}
};

} // end namespace

#endif