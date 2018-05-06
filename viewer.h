#ifndef _VIEWER_H_
#define _VIEWER_H_

#include "common.h"
#include "view.h"
#include "events.h"

namespace Renzoku {

class Viewer {
public:
    Viewer(int height, int width);
    ~Viewer();

public:
    /**
     * Start viewer. The first view is displayed by default if no previous views are 
     * explicitly tagged to display (by method \ref show(int)). 
     */
    virtual void run() = 0;

    /**
     * Schedule to redraw the views.
     */ 
    virtual void redisplay() = 0;

    /**
     * Add a view into the list of views of the scene. 
     * The View::init() function will be called by default.
     */
    void add(View *view);

    /**
     * Tag a specific view for display.
     */
    void show(int index); 
    
    /**
     * Issue requests to all child views to render from scratch as scene information
     * changes, e.g., camera parameters are updated. 
     */
    void reset_all_views();

    inline void set_scene(Scene *scene);

    inline int get_height() const;
    inline int get_width() const;

    inline View *get_cur_view() const;

protected:
    typedef vector<View *> Views;
    Views views;                        /// a list of views of the scene
    View *cur_view;                     /// the view that is being displayed
    int cur_view_index;                 /// index of the view being displayed (for view switch)
    int height, width;
    Scene *scene;

public:
    /**
     * GLUT-like functions
     */
    virtual void init(int argc, char **argv);
    virtual void reshape(int width, int height);
    virtual void timer(int val);    
    virtual void display();

    /**
     * Mouse down.
     */
    virtual void mouse(MouseButton button, MouseState state, int x, int y);

    /** 
     * Mouse move (with some buttons down).
     */
    virtual void motion(int x, int y);
    virtual void keyboard(unsigned char key, int x, int y);
};

inline int Viewer::get_width() const {
    return width;
}

inline int Viewer::get_height() const {
    return height;
}

inline void Viewer::set_scene(Scene *scene) {
    this->scene = scene;
}

inline View *Viewer::get_cur_view() const {
    return cur_view;
}

} // end namespace Renzoku

#endif

