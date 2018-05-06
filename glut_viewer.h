#ifndef _GLUT_VIEWER_H_
#define _GLUT_VIEWER_H_

#include "viewer.h"

namespace Renzoku {

class GlutViewer : public Viewer {
public:
    GlutViewer(int height, int width);    
    
    ~GlutViewer();

public:    
    void run();    
    void redisplay();

    /**
     * GLUT callback functions
     */
    void init(int argc, char **argv);
    void reshape(int width, int height);
    void timer(int val);    
    void display();
    void mouse(int button, int state, int x, int y);
    void motion(int x, int y);
    void keyboard(unsigned char key, int x, int y);
};

class GlutWrapper {
public:
    GlutWrapper();
    ~GlutWrapper();

public:
    void set_viewer(GlutViewer *viewer);
    
public:
    static void reshape(int width, int height);
    static void timer(int val);    
    static void display();
    static void mouse(int button, int state, int x, int y);
    static void motion(int x, int y);
    static void keyboard(unsigned char key, int x, int y);

    static MouseButton glut_to_mouse_button(int button);
    static MouseState  glut_to_mouse_state(int state);

private:
    static GlutViewer *viewer;
};

} // end namespace

#endif