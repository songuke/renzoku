#ifndef _VIEW_H_
#define _VIEW_H_

#include "common.h"
#include "events.h"

namespace Renzoku {

class Viewer; // forward declaration

/**
 * \brief View is an abstract class that defines the interface for scene views. 
 * A scene can be observed by different views, for example:
 * NullView (simply coordinate the renderer and forward the output to storage if any), 
 * GLView (MeshView, or ProgressiveView, which outputs
 * the scene to OpenGL framebuffer).
 *
 * View can adjust the rendering content based on simple user inputs from keyboard and mouse.
 */
class View {
public:
	View(Scene *scene) : scene(scene) {}

    /**
     * Called only once when the view is attached to the viewer.
     */
    virtual void init() {}

    /**
     * Clear all existing computed results and perform rendering from scratch.
     */
    virtual void reset() = 0;

    virtual void reshape(int width, int height) {}

	virtual void display() = 0;

    virtual void mouse(MouseButton button, MouseState state, int x, int y) {}
    virtual void motion(int x, int y) {}
    virtual void keyboard(unsigned char key, int x, int y) {}

    /**
     * Event fired when a view is switched to.
     */
    virtual void on_show() {}
   
    /**
     * Event fired when a view is deactivated and another view will be shown.
     */
    virtual void on_hide() {}

    /**
     * Event fired when the application closes.
     */
    virtual void on_close() {}

    void set_viewer(Viewer *viewer) { this->viewer = viewer; }

protected:
    Scene *scene;

    Viewer *viewer;
};

} // end namespace 

#endif
