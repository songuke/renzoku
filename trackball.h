// 15/05/2013: Trackball modified for proper implementation of zoom in projection matrix
//             and panning to multiply before gluLookat is called by Son Hua.
//             This fits more closely to the physical camera model.
// 10/29/2011: Trackball class removed for ANSI C compatibility by Son Hua. 
//  6/18/2000: Rotation animation disabled by Kok-Lim Low.
//  6/17/2000: Modified by Kok-Lim Low to include panning and zooming.

#ifndef _TRACKBALL_H_
#define _TRACKBALL_H_

#include "events.h"

/* 
 *  Simple trackball-like motion adapted (ripped off) from projtex.c
 *  (written by David Yu and David Blythe).  See the SIGGRAPH '96
 *  Advanced OpenGL course notes.
 *
 *  This trackball works with the following transformation order (read from bottom to top
 *  as OpenGL convention): 
 *  -- Load identity for modelview matrix
 *  -- Apply panning
 *  -- Translate world back to default
 *  -- Apply rotation and zoom
 *  -- Translate world to center at the centroid of the scene
 *  -- Draw objects
 *
 *  Usage:
 *
 *  o  create a TrackBall object, say called tb
 *  o  call tb_reset() to reset to the initial transformation
 *  o  call tb_reshape() from the reshape callback
 *  o  call tb_apply_rotation_zoom() and tb_apply_panning() to apply the trackball transformation. See the order above.
 *  o  call tb_motion() from the motion callback
 *  o  call tb_mouse() from the mouse callback
 *  o  note that tb_reshape() must be called to update 
 *     width and height of the viewport for rotation to work. 
 *
 *  Typical setup:
 
	#include "glext.h"
    #include "trackball.h"


	TrackBall tb(GLUT_LEFT_BUTTON, GLUT_MIDDLE_BUTTON, GLUT_RIGHT_BUTTON);


	void display(void)
    {
		. . .
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

        tb_apply_panning();

        gluLookat(...) or glTranslatef(0.0, 0.0, -viewer_distance_from_object_center);
        
		glPushMatrix();
            glTranslatef(obj_center_x, obj_center_y, obj_center_z);
			tb_apply_rotation_zoom();
			glTranslatef(-obj_center_x, -obj_center_y, -obj_center_z);

			. . . draw the scene . . .

		glPopMatrix();
		. . .
    } 

    void reshape(int width, int height)
    {
		tb_reshape(width, height);
    }

    void mouse(int button, int state, int x, int y)
    {
		tb_mouse(button, state, x, y);
		glutPostRedisplay();
    }


    void motion(int x, int y)
    {
		tb_motion(x, y);
		glutPostRedisplay();
    } 


    int main(int argc, char** argv)
    {		
		glutDisplayFunc(display);
		glutReshapeFunc(reshape);
		glutMouseFunc(mouse);
		glutMotionFunc(motion);		
    }
 *
 * */

void tb_init(void);
void tb_init_buttons(MouseButton rot_button, MouseButton pan_button, MouseButton zoom_button);	
void tb_reset(void);

/**
 * Panning moves current eye and lookat along the right and up vector (if we think in world coordinates). 
 * In camera coordinates, it is simply a translation along the X and Y axes of the camera coordinate system.
 * 
 * This function implements panning in eye space. 
 * Since OpenGL performs matrix post-multiplication, this function needs to be called before gluLookat() or 
 * eye transformation.
 */
void tb_apply_panning(void);
void tb_apply_rotation_zoom(void);

void tb_reshape(int width, int height);
void tb_mouse(MouseButton button, MouseState state, int x, int y);
void tb_motion(int x, int y);
void tb_load(const char *file);
void tb_save(const char *file);

/**
 * Return zoom and pan parameters for physical camera model update if any. 
 */
GLdouble tb_get_zoom(void);
void tb_get_pan(GLdouble pan[2]);

#endif
