#include "viewer.h"
#include "log.h"

namespace Renzoku {

Viewer::Viewer(int height, int width) : height(height), width(width), cur_view(NULL), scene(NULL) {

}

Viewer::~Viewer() {
    
}

void Viewer::add(View *view) {    
    views.push_back(view);
    
    view->set_viewer(this);
    view->init();
}

void Viewer::display() {
    cur_view->display();
}

void Viewer::reshape(int w, int h) {
    this->width = w;
    this->height = h;
    
    for (int i = 0; i < views.size(); ++i)
        views[i]->reshape(w, h);
}

void Viewer::mouse(MouseButton button, MouseState state, int x, int y) {
    cur_view->mouse(button, state, x, y);
}

void Viewer::motion(int x, int y) {
    cur_view->motion(x, y);
}

void Viewer::keyboard(unsigned char key, int x, int y) {  
    switch (key) {
        case ' ':
            cur_view->keyboard(key, x, y);

            // switch to next view            
            this->show((cur_view_index + 1) % views.size());

            break;

        default:
            cur_view->keyboard(key, x, y);
            break;
    }  
}

void Viewer::timer(int val) {
}

void Viewer::init(int argc, char **argv) {
}

void Viewer::show(int index) {
    if (cur_view)
        cur_view->on_hide();

    cur_view = views[index];
    cur_view_index = index;

    cur_view->on_show();
}

void Viewer::run() {
}

void Viewer::reset_all_views() {
    for (int i = 0; i < views.size(); ++i)
        views[i]->reset();
}

} // end namespace Renzoku
