#ifndef _MESH_VIEW_H_
#define _MESH_VIEW_H_

#include "gl_view.h"
#include "boundingbox.h"
#include "camera.h"
#include "mutable_path.h"

#include <boost/thread.hpp>

namespace Renzoku {

struct MeshDisplay {
    enum Mode {
        WIREFRAME,
        FILL,
        NONE,
        NUM_DISPLAY_MODES
    };
};

class MeshView : public GLView {
public:
    MeshView(Scene *scene);
    MeshView(Scene *scene, int height, int width);
    ~MeshView();

    void init();
    void reset();

    void reshape(int width, int height);
    void display();
    void keyboard(unsigned char key, int x, int y);
    void mouse(MouseButton button, MouseState state, int x, int y);
    void motion(int x, int y);

    void on_show();
    void on_hide();

    /**
     * For debugging.
     */
    void set_bounding_boxes(BoundingBoxes &boxes);
    void set_path(const MutablePath &path);
    void set_draw_paths(bool draw);

protected:
    void draw_axes(float length) const;
    void draw_bounding_box(BoundingBox *box, Rgb color = DefaultRgb::grey) const;
    void draw_canonical_bounding_box() const;
    void draw_path(const MutablePath &path, Rgb color = DefaultRgb::white) const;

    void update_camera_from_trackball();

    /** 
     * Read scene camera parameters and apply to internal camera parameters. 
     */
    void update_camera_from_scene();

    /**
     * Set OpenGL camera with current eye, lookat, up, and field of view.
     */
    void reset_camera();
    
protected:
    int height, width;                  /// derive from the size of the viewer
    Vec3 eye, lookat, up;
    Float fov, near_plane, far_plane;

    /*
    GLuint vbo;
    GLuint ibo;
    GLuint vbo_offset;
    GLuint ibo_offset; 
    const GLuint vbo_max_size = 8*2048*2048; // 32 MB (abt. 800k vertices)
    const GLuint ibo_max_size = 8*2048*2048; 

    GLuint fbo;
    GLuint fbo_tex;
    GLuint fbo_depth;*/

    MeshDisplay::Mode display_mode;

	BoundingBoxes boxes;                 /// spatial partition bounding boxes, e.g., from an Octree. For debugging purpose. 
	bool is_draw_spatial_bounding_boxes;

    vector<MutablePath> paths;
    bool is_draw_paths;
    boost::mutex mtx_paths;
};

} // end namespace

#endif
