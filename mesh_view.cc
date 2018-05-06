#include "noglut.h"

#include "mesh_view.h"
#include "scene.h"
#include "aggregate.h"
#include "boundingbox.h"

#include "viewer.h"

#include "triangle.h"
#include "quad.h"
#include "sphere.h"

#include "area_light.h"
#include "point_light.h"

extern "C" {
    #include "trackball.h"
}

namespace Renzoku {

extern "C" {
void APIENTRY
glutWireSphere(GLdouble radius, GLint slices, GLint stacks);
void APIENTRY
glutSolidSphere(GLdouble radius, GLint slices, GLint stacks);
}

MeshView::MeshView(Scene *scene) : GLView(scene), height(512), width(512) {
    display_mode = MeshDisplay::WIREFRAME;

	is_draw_spatial_bounding_boxes = false;
    is_draw_paths = false;
    scene->set_mesh_view(this);
}

MeshView::MeshView(Scene *scene, int height, int width) : GLView(scene), height(height), width(width) {
    display_mode = MeshDisplay::WIREFRAME;

	is_draw_spatial_bounding_boxes = false;
    is_draw_paths = false;
    scene->set_mesh_view(this);
}

MeshView::~MeshView() {

}

void MeshView::init() {
    BoundingBox box = scene->get_aggregate()->get_bounding_box();
    Vec3 centroid = box.centroid();
    Vec3 size = box.size();
    
    eye = centroid + 2 * (box.v_max - centroid);      // simulate a view from (1, 1, 1) to (0, 0, 0)
    lookat = centroid;
    up = Vec3(0, 1, 0);
    
    fov = 60.0;

    Camera *camera = scene->get_camera();
    if (camera->get_near_plane() > 0) {
        near_plane = camera->get_near_plane();
    } else {
        Float min_near_plane = size.min_component() * 0.05f;     // for good depth resolution, the point light won't see anything nearer than 5% of the minimum bounding box size.
        near_plane = std::max(min_near_plane, 0.05f * size.max_component());    
    }
       
    if (camera->get_far_plane() > 0) {
        far_plane = camera->get_far_plane();
    } else {
        far_plane = near_plane + 10 * size.max_component();
    }

    this->reset_camera();
        
    tb_init_buttons(MOUSE_BUTTON_LEFT, MOUSE_BUTTON_MIDDLE, MOUSE_BUTTON_RIGHT);
}

void MeshView::reset_camera() {
    // set camera with current eye and lookat.
    // Trackball status is preserved as before hiding.
    glViewport(0, 0, width, height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(fov, width * 1.0 / height, near_plane, far_plane);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(eye.x(), eye.y(), eye.z(),
              lookat.x(), lookat.y(), lookat.z(),
              up.x(), up.y(), up.z());
}

void MeshView::on_show() {
    update_camera_from_scene();
    //reset_camera();
}

void MeshView::on_hide() {
}

void MeshView::reshape(int width, int height) {
    this->width = width;
    this->height = height;

    tb_reshape(width, height);
    glViewport(0, 0, width, height);
}

/**
 * FIXME: this violates LSP design principle. 
 */
static void draw_shape(Shape *s, const Rgb &color) {
    glColor3f(color.red(), color.green(), color.blue());
        
    Quad *q = dynamic_cast<Quad *>(s); // can have overhead due to run-time type check
    if (q) {
        glBegin(GL_QUADS);
        glVertex3f(q->p0.x(), q->p0.y(), q->p0.z());
        glVertex3f(q->p1.x(), q->p1.y(), q->p1.z());
        glVertex3f(q->p2.x(), q->p2.y(), q->p2.z());
        glVertex3f(q->p3.x(), q->p3.y(), q->p3.z());
        glEnd();
    } else {
        Triangle *t = dynamic_cast<Triangle *>(s);
        if (t) {
            glBegin(GL_TRIANGLES);
            glVertex3f(t->p0.x(), t->p0.y(), t->p0.z());
            glVertex3f(t->p1.x(), t->p1.y(), t->p1.z());
            glVertex3f(t->p2.x(), t->p2.y(), t->p2.z());
            glEnd();
        } else {
            Sphere *e = dynamic_cast<Sphere *>(s);
            if (e) {
                glPushMatrix();
                glTranslatef(e->center.x(), e->center.y(), e->center.z());
                glutSolidSphere(e->rad, 16, 16);
                glPopMatrix();
            } else {
                // unknown shape
            }
        }
    }
}
    
void MeshView::display() {
    BoundingBox box = scene->get_aggregate()->get_bounding_box();
    Vec3 centroid = box.centroid();
    Vec3 size = box.size();

    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);        
        
    // ----- Trackball usage -------------------------------------------------------------------------
    // Trackball allows rotation around the centroid of the scene. Therefore,     
    // we translate to the centroid, perform rotation (and zoom), and move back to original world coordinates. 
    // Then we can apply eye view transformation, and finally panning. 
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    tb_apply_panning();

    gluLookAt(eye.x(), eye.y(), eye.z(),
              lookat.x(), lookat.y(), lookat.z(),
              up.x(), up.y(), up.z());              
    
    glTranslatef(centroid.x(), centroid.y(), centroid.z());
    tb_apply_rotation_zoom();
    glTranslatef(-centroid.x(), -centroid.y(), -centroid.z());
    // -----------------------------------------------------------------------------------------------

    // draw axes at the center of the bounding box
    glPushMatrix();
    glTranslatef(centroid.x(), centroid.y(), centroid.z());
    draw_axes(box.size().max_component());
    glPopMatrix();

    draw_bounding_box(&box);
	
	// draw boxes from spatial partition data structure if requested
	if (is_draw_spatial_bounding_boxes) {
		for (int i = 0; i < boxes.size(); ++i)
			draw_bounding_box(&boxes[i], DefaultRgb::green);
	}		

    if (is_draw_paths) {
        if (paths.size() > 0) {
            MutablePath path(scene);
            mtx_paths.lock();
            path = paths[0];
            mtx_paths.unlock();
            draw_path(path, DefaultRgb::white);
        }
    }

    if (display_mode == MeshDisplay::NONE) return;

    glEnable(GL_DEPTH_TEST);
    //glDisable(GL_CULL_FACE);
    switch (display_mode) {
    case MeshDisplay::WIREFRAME:                
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        break;
    case MeshDisplay::FILL:
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        break;
    }

    Surfaces& surfaces = *scene->get_surfaces();
    for (int i = 0; i < surfaces.size(); ++i) {
        Rgb color;
        if (surfaces[i].is_light()) {
            color = surfaces[i].get_area_light()->power();
        } else if (surfaces[i].is_env_light()) {
            color = DefaultRgb::white;
        } else {
            color = surfaces[i].get_material()->get_representative_color();
        }

        Shape *s = surfaces[i].get_shape();
        
        draw_shape(s, color);
    }

    Lights& lights = *scene->get_lights();
    for (int i = 0; i < lights.size(); ++i) {
        if (lights[i]->get_light_type() == Light::AREA_LIGHT) continue;

        if (lights[i]->get_light_type() == Light::POINT_LIGHT) {
            PointLight *pl = (PointLight *)lights[i];
            Rgb color = pl->power();
            Vec3 pos = pl->org();
            glColor3f(color.red(), color.green(), color.blue());
            glPushMatrix();
            glTranslatef(pos.x(), pos.y(), pos.z());
            glutSolidSphere(8, 16, 16);
            glPopMatrix();
        }
    }
}

/**
 * Calculate eye and lookat after applying trackball transformation
 */
void MeshView::update_camera_from_trackball() {
    GLdouble m[4*4]; 
    BoundingBox box = scene->get_aggregate()->get_bounding_box();
    Vec3 centroid = box.centroid();
    Vec3 size = box.size();

    // query the camera matrix (without panning) to solve for eye and lookat, and apply panning later

    // simulate trackball transform to obtain the view matrix (without panning)
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();

    glLoadIdentity();
    gluLookAt(eye.x(), eye.y(), eye.z(),
              lookat.x(), lookat.y(), lookat.z(),
              up.x(), up.y(), up.z());  
    glTranslatef(centroid.x(), centroid.y(), centroid.z());
    tb_apply_rotation_zoom();
    glTranslatef(-centroid.x(), -centroid.y(), -centroid.z());
        
    glGetDoublev(GL_MODELVIEW_MATRIX, m);
    glPopMatrix();

    // zoom is implemented in projection matrix, so UVN are unit by default.
    Vec3 u(m[0], m[4], m[8]); //u.make_unit_vector();   
    Vec3 v(m[1], m[5], m[9]); //v.make_unit_vector();
    Vec3 n(m[2], m[6], m[10]); //n.make_unit_vector();
    Vec3 dot_eye(-m[12], -m[13], -m[14]);
    eye = Vec3(u.x() * dot_eye.x() + v.x() * dot_eye.y() + n.x() * dot_eye.z(),
               u.y() * dot_eye.x() + v.y() * dot_eye.y() + n.y() * dot_eye.z(),
               u.z() * dot_eye.x() + v.z() * dot_eye.y() + n.z() * dot_eye.z());

    lookat = eye - n * size.max_component();
    up = v; 

    // pan
    GLdouble pan[2];
    tb_get_pan(pan);
    eye    += -pan[0] * u + -pan[1] * v;
    lookat += -pan[0] * u + -pan[1] * v;

    GLdouble zoom = tb_get_zoom();
    fov = atan(tan(fov / 2 * M_PI / 180.0) / zoom) * 180.0 / M_PI * 2;

    Camera *camera = scene->get_camera();
    Size2 film_size = camera->get_film_size();
    Float focal = Camera::focal_length_from_fov(fov, film_size.height);
    cout << "Eye    : " << eye << endl;
    cout << "Lookat : " << lookat << endl;
    cout << "Up     : " << up << endl;
    cout << "Near   : " << near_plane << endl;
    cout << "Far    : " << far_plane << endl;
    cout << "Zoom   : " << zoom << endl;
    cout << "FOV    : " << fov << endl;
    cout << "Focal  : " << focal << endl;
    cout << "Film (w x h)   : " << film_size.width << " " << film_size.height << endl;

    // after applying trackball transform, reset trackball to identity
    reset_camera();
    tb_reset();
}

void MeshView::update_camera_from_scene() {
    Camera *camera = scene->get_camera();
    eye = camera->get_eye();
    lookat = camera->get_lookat();
    up = camera->get_up();
    fov = camera->get_vertical_fov();

    // sync with OpenGL camera
    reset_camera();
}

void MeshView::keyboard(unsigned char key, int x, int y) {
    switch (key) {
        case 't':   // test if the trackball and current FOV calculation is consistent
            update_camera_from_trackball();            
            viewer->redisplay();
            break;
            
        case ' ':            // event from Viewer to change view
        case 'c':            // just update, no change view
        {
            // apply current camera settings to all views
            // by changing the camera stored in Scene.
            Camera *c = scene->get_camera()->clone();
            update_camera_from_trackball();
            c->set_perspective(eye, lookat, up, fov);

            scene->set_camera(c);
                        
            viewer->reset_all_views();
            break;
        }

        case 'w':
        {
            update_camera_from_trackball();

            // move the eye forward 10% of the bounding box size
            Float size = scene->get_aggregate()->get_bounding_box().size().max_component();
            eye += unit_vector(lookat - eye) * size * 0.1f; 

            viewer->redisplay();
            break;
        }

        case 's':
        {
            update_camera_from_trackball();

            // move the eye backward 10% of the bounding box size
            Float size = scene->get_aggregate()->get_bounding_box().size().max_component();
            eye -= unit_vector(lookat - eye) * size * 0.1; 

            viewer->redisplay();
            break;
        }

        case 'm':
            display_mode = (MeshDisplay::Mode)((display_mode + 1) % MeshDisplay::NUM_DISPLAY_MODES);
            viewer->redisplay();
            break;

    	case 'b':
			is_draw_spatial_bounding_boxes ^= 1;
			viewer->redisplay();
			break;

        // TODO: add a hot key to change view to see the whole scene (original trackball view).
        case 'd':
            // default bounding view (for cases where camera of the scene is wrongly set)
            BoundingBox box = scene->get_aggregate()->get_bounding_box();
            Vec3 centroid = box.centroid();
            Vec3 size = box.size();
            eye = centroid + 2 * (box.v_max - centroid);      // simulate a view from (1, 1, 1) to (0, 0, 0)
            lookat = centroid;
            up = Vec3(0, 1, 0);    
            fov = 60.0;
            scene->get_camera()->set_perspective(eye, lookat, up, fov);
            viewer->redisplay();
            break;
    }
}

void MeshView::mouse(MouseButton button, MouseState state, int x, int y) {        
    tb_mouse(button, state, x, y);
    viewer->redisplay();
}
    
void MeshView::motion(int x, int y) {
    tb_motion(x, y);
    viewer->redisplay();
}

void MeshView::reset() {
    
}

void MeshView::draw_axes(float length) const {
	glPushAttrib(GL_ALL_ATTRIB_BITS);
    glDisable(GL_LIGHTING);
    glLineWidth(3.0);
    glBegin(GL_LINES);
    // x-axis.
    glColor3f(1.f, 0.f, 0.f);
    glVertex3f(0.f, 0.f, 0.f );
    glVertex3f(length, 0.f, 0.f);
    // y-axis.
    glColor3f(0.f, 1.f, 0.f);
    glVertex3f(0.f, 0.f, 0.f);
    glVertex3f(0.f, length, 0.f);
    // z-axis.
    glColor3f(0.f, 0.f, 1.f);
    glVertex3f(0.f, 0.f, 0.f);
    glVertex3f(0.f, 0.f, length);
    glEnd();
    glPopAttrib();
}

/**
 * Draw a cube by lines from (0, 0, 0) to (1, 1, 1)
 */
void MeshView::draw_canonical_bounding_box() const {
    glBegin(GL_LINES);
	glVertex3f(0, 0, 0);
	glVertex3f(1, 0, 0);
	glVertex3f(0, 0, 0);
	glVertex3f(0, 1, 0);
	glVertex3f(0, 0, 0);
	glVertex3f(0, 0, 1);

	glVertex3f(1, 1, 0);
	glVertex3f(0, 1, 0);
	glVertex3f(1, 1, 0);
	glVertex3f(1, 0, 0);
	glVertex3f(1, 1, 0);
	glVertex3f(1, 1, 1);

	glVertex3f(1, 0, 1);
	glVertex3f(0, 0, 1);
	glVertex3f(1, 0, 1);
	glVertex3f(1, 1, 1);
	glVertex3f(1, 0, 1);
	glVertex3f(1, 0, 0);

	glVertex3f(0, 1, 1);
	glVertex3f(1, 1, 1);
	glVertex3f(0, 1, 1);
	glVertex3f(0, 0, 1);
	glVertex3f(0, 1, 1);
	glVertex3f(0, 1, 0);
    glEnd();
}

void MeshView::draw_bounding_box(BoundingBox *box, Rgb color) const {
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glDisable(GL_LIGHTING);
    glLineWidth(1.0);
	glColor4f(color.r, color.g, color.b, 1.0f);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();

    Vec3 size = box->size();
    glTranslatef(box->v_min.x(), box->v_min.y(), box->v_min.z());
    glScalef(size.x(), size.y(), size.z());    

    draw_canonical_bounding_box();
    
    glPopMatrix();
    glPopAttrib();
}

void MeshView::set_bounding_boxes(BoundingBoxes &b) {
    boxes.clear();
    boxes.assign(b.begin(), b.end());
}

void MeshView::draw_path(const MutablePath &path, Rgb color) const {
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glDisable(GL_LIGHTING);
    glLineWidth(1.0);
    glColor4f(color.r, color.g, color.b, 1.0f);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();

    glBegin(GL_LINE_STRIP);
    for (int k = 0; k < path.num_nodes; ++k) {
        Vec3 p = path.nodes[k].dg.p;
        glVertex3f(p.x(), p.y(), p.z());
    }
    glEnd();

    glPopMatrix();
    glPopAttrib();
}

void MeshView::set_path(const MutablePath &path) {
    mtx_paths.lock();

    paths.clear();
    paths.push_back(path);

    mtx_paths.unlock();
}

void MeshView::set_draw_paths(bool draw) {
    is_draw_paths = draw;
}

// ----------------------------------------------------------------------------
// Some shape drawing functions from GLUT 3.7 source code.
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// glut_shapes.c
// ----------------------------------------------------------------------------
static GLUquadricObj *quadObj;

#define QUAD_OBJ_INIT() { if(!quadObj) initQuadObj(); }

static void
initQuadObj(void)
{
  quadObj = gluNewQuadric();
  //if (!quadObj)
  //  __glutFatalError("out of memory.");
}

/* CENTRY */
void APIENTRY
glutWireSphere(GLdouble radius, GLint slices, GLint stacks)
{
  QUAD_OBJ_INIT();
  gluQuadricDrawStyle(quadObj, GLU_LINE);
  gluQuadricNormals(quadObj, GLU_SMOOTH);
  /* If we ever changed/used the texture or orientation state
     of quadObj, we'd need to change it to the defaults here
     with gluQuadricTexture and/or gluQuadricOrientation. */
  gluSphere(quadObj, radius, slices, stacks);
}

void APIENTRY
glutSolidSphere(GLdouble radius, GLint slices, GLint stacks)
{
  QUAD_OBJ_INIT();
  gluQuadricDrawStyle(quadObj, GLU_FILL);
  gluQuadricNormals(quadObj, GLU_SMOOTH);
  /* If we ever changed/used the texture or orientation state
     of quadObj, we'd need to change it to the defaults here
     with gluQuadricTexture and/or gluQuadricOrientation. */
  gluSphere(quadObj, radius, slices, stacks);
}

} // end namespace
