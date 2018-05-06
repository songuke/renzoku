#include <GL/glew.h>
#include <GL/freeglut.h> // use FreeGLUT instead of GLUT for on-closing event

#include "glut_viewer.h"
#include "scene.h"
#include "stats.h"
#include "log.h"

#include <sstream>
using namespace std;

namespace Renzoku {

// define the static viewer variable of GlutWrapper here. Otherwise undefined reference linking.
GlutViewer *GlutWrapper::viewer = NULL;
static GlutWrapper gw;

GlutWrapper::GlutWrapper() {
}

GlutWrapper::~GlutWrapper() {
}

void GlutWrapper::set_viewer(GlutViewer *_viewer) {
    viewer = _viewer;
}

void GlutWrapper::reshape(int width, int height) {
    viewer->reshape(width, height);
}

void GlutWrapper::timer(int val) {
    viewer->timer(val);
}

void GlutWrapper::display() {
    viewer->display();
}

void GlutWrapper::mouse(int button, int state, int x, int y) {
    viewer->mouse(button, state, x, y);
}

void GlutWrapper::motion(int x, int y) {
    viewer->motion(x, y);
}

void GlutWrapper::keyboard(unsigned char key, int x, int y) {
    viewer->keyboard(key, x, y);
}

MouseButton GlutWrapper::glut_to_mouse_button(int button) {
    switch (button) {
    case GLUT_LEFT_BUTTON:
        return MOUSE_BUTTON_LEFT;
    case GLUT_RIGHT_BUTTON:
        return MOUSE_BUTTON_RIGHT;
    case GLUT_MIDDLE_BUTTON:
        return MOUSE_BUTTON_MIDDLE;
    }
    return MOUSE_BUTTON_LEFT;
}

MouseState GlutWrapper::glut_to_mouse_state(int state) {
    switch (state) {
    case GLUT_DOWN:
        return MOUSE_STATE_DOWN;
    case GLUT_UP:
        return MOUSE_STATE_UP;
    case GLUT_DOUBLE:
        return MOUSE_STATE_DOUBLE;
    }
    return MOUSE_STATE_UP;
}

GlutViewer::GlutViewer(int height, int width) : Viewer(height, width) {
    gw.set_viewer(this);
}

GlutViewer::~GlutViewer() {
    
}

void GlutViewer::display() {
    Viewer::display();

    ostringstream oss;
    oss << "FPS: " << scene->get_stats_counter()->fps();
    glutSetWindowTitle(oss.str().c_str());

    glutSwapBuffers();
}

void GlutViewer::reshape(int w, int h) {
    Viewer::reshape(w, h);
}

void GlutViewer::mouse(int button, int state, int x, int y) {
    MouseButton m_button = gw.glut_to_mouse_button(button);
    MouseState m_state   = gw.glut_to_mouse_state(state);

    Viewer::mouse(m_button, m_state, x, y);
}

void GlutViewer::motion(int x, int y) {
    Viewer::motion(x, y);
}

void GlutViewer::keyboard(unsigned char key, int x, int y) {
    Viewer::keyboard(key, x, y);  
}

void GlutViewer::timer(int val) {
    glutPostRedisplay();
    glutTimerFunc(1, GlutWrapper::timer, 0);
}

void GlutViewer::init(int argc, char **argv) {
    //glutInitContextVersion (3, 2);
    //glutInitContextFlags (GLUT_FORWARD_COMPATIBLE);
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);    
    glutInitWindowSize(width, height); 
    glutCreateWindow("Renzoku");

    glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS); 

    glutReshapeFunc(GlutWrapper::reshape);
    glutDisplayFunc(GlutWrapper::display);
    glutTimerFunc(1, GlutWrapper::timer, 0);  
    glutMouseFunc(GlutWrapper::mouse);
    glutMotionFunc(GlutWrapper::motion);
    glutKeyboardFunc(GlutWrapper::keyboard);
        
    /*
    GLenum err = glewInit();
    if (GLEW_OK != err) {
        string str((char *)glewGetErrorString(err));
        Log::critical(str, ExitCode::EXIT_CODE_LIBRARY_FAILED_TO_LOAD);
    }
    Log::info() << "GLEW: " << glewGetString(GLEW_VERSION) << endl;
    */
}

void GlutViewer::run() {
    if (cur_view == NULL) {        
        show(0);
    }

    glutMainLoop();    

    for (int i = 0; i < views.size(); ++i)
        views[i]->on_close();    
}

void GlutViewer::redisplay() {
    glutPostRedisplay();
}

} // end namespace