//  6/18/2000: Rotation animation disabled by Kok-Lim Low.
//  6/17/2000: Modified by Kok-Lim Low to include panning and zooming.

/*
 *  Simple trackball-like motion adapted (ripped off) from projtex.c
 *  (written by David Yu and David Blythe).  See the SIGGRAPH '96
 *  Advanced OpenGL course notes.
 */

#include "noglut.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "trackball.h"

static struct {
	GLdouble tb_angle;
	GLdouble tb_axis[3];
	GLdouble tb_transform[4*4];

	GLuint tb_width;
	GLuint tb_height;

	GLdouble tb_pan_x;
	GLdouble tb_pan_y;    

	GLdouble tb_zoom;           /// one-time zoom factor which is accumulated to the projection matrix directly. This 
                                /// Therefore, this value is reset after it is applied.
	GLdouble tb_zoom_inc;
    GLdouble tb_zoom_all;       /// the accumulated zoom factor since the original focal length, for retrieving new field of view 
                                /// and focal length.

	GLint tb_rot_button;
	GLint tb_pan_button;
	GLint tb_zoom_button;

	GLint tb_mouse_button;
	GLint tb_mouse_x;
	GLint tb_mouse_y;

	GLdouble tb_model_mat[4*4];
	GLdouble tb_proj_mat[4*4];
	GLint tb_viewport[4];
} tb;

static void _tbPointToVector(int x, int y, int width, int height, GLdouble v[3])
{
	GLdouble d, a;

	// project x, y onto a hemi-sphere centered within width, height.
	v[0] = (2.0 * x - width) / width;
	v[1] = (height - 2.0 * y) / height;
	d = sqrt(v[0] * v[0] + v[1] * v[1]);
	v[2] = cos((3.14159265 / 2.0) * ((d < 1.0) ? d : 1.0));
	a = 1.0 / sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	v[0] *= a;
	v[1] *= a;
	v[2] *= a;
}

static void tb_capture_transform(void)
{
	glGetDoublev(GL_MODELVIEW_MATRIX, tb.tb_model_mat);
	glGetDoublev(GL_PROJECTION_MATRIX, tb.tb_proj_mat);
	glGetIntegerv(GL_VIEWPORT, tb.tb_viewport);
}

void tb_init() 
{	
	tb.tb_rot_button = MOUSE_BUTTON_LEFT;
	tb.tb_pan_button = MOUSE_BUTTON_MIDDLE;
	tb.tb_zoom_button = MOUSE_BUTTON_RIGHT;

	tb_reset();
    tb.tb_axis[0] = 1.0; tb.tb_axis[1] = 0.0; tb.tb_axis[2] = 0.0;
}

void tb_init_buttons(MouseButton rot_button, MouseButton pan_button, MouseButton zoom_button)
{
	if (rot_button == pan_button || rot_button == zoom_button || pan_button == zoom_button)
	{
		tb.tb_rot_button = MOUSE_BUTTON_LEFT;
		tb.tb_pan_button = MOUSE_BUTTON_MIDDLE;
		tb.tb_zoom_button = MOUSE_BUTTON_RIGHT;
	}
	else
	{
		tb.tb_rot_button = rot_button;
		tb.tb_pan_button = pan_button;
		tb.tb_zoom_button = zoom_button;
	}

	tb_reset();
}



void tb_reset(void)
{
	tb.tb_angle = 0.0;
	//tb.tb_axis[0] = 1.0; tb.tb_axis[1] = 0.0; tb.tb_axis[2] = 0.0;

    tb.tb_pan_x = 0.0;
	tb.tb_pan_y = 0.0;    
    
    tb.tb_zoom = 1.0;
	tb.tb_zoom_inc = 0.005;
    tb.tb_zoom_all = 1.0;

	// put the identity in the trackball transform

	tb.tb_transform[0] = 1.0; tb.tb_transform[4] = 0.0; tb.tb_transform[8] = 0.0; tb.tb_transform[12] = 0.0;
	tb.tb_transform[1] = 0.0; tb.tb_transform[5] = 1.0; tb.tb_transform[9] = 0.0; tb.tb_transform[13] = 0.0;
	tb.tb_transform[2] = 0.0; tb.tb_transform[6] = 0.0; tb.tb_transform[10] = 1.0; tb.tb_transform[14] = 0.0;
	tb.tb_transform[3] = 0.0; tb.tb_transform[7] = 0.0; tb.tb_transform[11] = 0.0; tb.tb_transform[15] = 1.0;
}


void tb_apply_rotation_zoom(void)
{
    // prepare rotation matrix
	glPushMatrix();
	glLoadIdentity();
	glRotated(tb.tb_angle, tb.tb_axis[0], tb.tb_axis[1], tb.tb_axis[2]);
	glMultMatrixd(tb.tb_transform);
	glGetDoublev(GL_MODELVIEW_MATRIX, tb.tb_transform);
	glPopMatrix();

    // Originally, this is the original hack that performs zoom in camera space. 
    // To make it "physically correct", zoom should be performed with "focal length", or "fov". 
    // which means it should be performed together with projection matrix.
    /*
	glTranslated(tb.tb_pan_x, tb.tb_pan_y, 0.0);   // perspective move the center of projection
	glScaled(tb.tb_zoom, tb.tb_zoom, tb.tb_zoom);  // scale the focal length
	glMultMatrixd(tb.tb_transform);
    */
        
    glMatrixMode(GL_PROJECTION); 	    
    glScaled(tb.tb_zoom, tb.tb_zoom, 1); // this is equivalent to calling gluPerspective again with the new FOV.
        
    glMatrixMode(GL_MODELVIEW);
    glMultMatrixd(tb.tb_transform);
        
    tb.tb_angle = 0.0;
	//tb.tb_axis[0] = 1.0; tb.tb_axis[1] = 0.0; tb.tb_axis[2] = 0.0;

    // need to reset zoom since projection matrix is set once only    
    tb.tb_zoom_all *= tb.tb_zoom;
    tb.tb_zoom = 1.0;
}

void tb_apply_panning(void) { 
    glMatrixMode(GL_MODELVIEW);   
    glTranslated(tb.tb_pan_x, tb.tb_pan_y, 0.0); 
}

void tb_reshape(int width, int height)
{
	tb.tb_width  = width;
	tb.tb_height = height;
}



void tb_mouse(MouseButton button, MouseState state, int x, int y)
{
	if (state == MOUSE_STATE_DOWN) tb.tb_mouse_button = button;
	tb.tb_mouse_x = x;
	tb.tb_mouse_y = y;
}



void tb_motion(int x, int y)
{
	GLdouble last_position[3], current_position[3], dx, dy, dz;
	GLdouble winx, winy, winz, tmp, old_pan_x, old_pan_y, new_pan_x, new_pan_y;

	tb.tb_angle = 0.0;
	tb.tb_axis[0] = 1.0; tb.tb_axis[1] = 0.0; tb.tb_axis[2] = 0.0;

	// rotating
	if (tb.tb_mouse_button == tb.tb_rot_button)
	{
		_tbPointToVector(tb.tb_mouse_x, tb.tb_mouse_y, tb.tb_width, tb.tb_height, last_position);
		_tbPointToVector(x, y, tb.tb_width, tb.tb_height, current_position);

		// calculate the angle to rotate by (directly proportional to the
		// length of the mouse movement
		dx = current_position[0] - last_position[0];
		dy = current_position[1] - last_position[1];
		dz = current_position[2] - last_position[2];
		tb.tb_angle = 90.0 * sqrt(dx * dx + dy * dy + dz * dz);

		// calculate the axis of rotation (cross product)
		tb.tb_axis[0] = last_position[1] * current_position[2] - 
			            last_position[2] * current_position[1];
		tb.tb_axis[1] = last_position[2] * current_position[0] - 
				        last_position[0] * current_position[2];
		tb.tb_axis[2] = last_position[0] * current_position[1] - 
				        last_position[1] * current_position[0];
	}


	// panning
	else if (tb.tb_mouse_button == tb.tb_pan_button)
	{        
        tb_capture_transform();

		gluProject(0.0, 0.0, 0.0, tb.tb_model_mat, tb.tb_proj_mat, tb.tb_viewport, &winx, &winy, &winz);
		gluUnProject(x, y, winz, tb.tb_model_mat, tb.tb_proj_mat, tb.tb_viewport, 
					 &new_pan_x, &new_pan_y, &tmp);
        
		gluUnProject(tb.tb_mouse_x, tb.tb_mouse_y, winz, tb.tb_model_mat, tb.tb_proj_mat, tb.tb_viewport, 
					  &old_pan_x, &old_pan_y, &tmp);
		
		tb.tb_pan_x += (new_pan_x - old_pan_x);
		tb.tb_pan_y -= (new_pan_y - old_pan_y);
	}


	// zooming or eye forwarding
	else if (tb.tb_mouse_button == tb.tb_zoom_button)
	{
        tb.tb_zoom += (tb.tb_mouse_y - y) * tb.tb_zoom_inc;
		if (tb.tb_zoom <= 0.0) tb.tb_zoom = tb.tb_zoom_inc;        
	}

	tb.tb_mouse_x = x;
	tb.tb_mouse_y = y;
}



void tb_save(const char *file)
{
    int i;
	FILE *fp = fopen(file, "w");	
    if (fp == NULL)
	{
		fprintf(fp, "TrackBall(): Cannot open file \"%s\" for writing.\n", file);
		return;
	}

	fprintf(fp, "%.10f\n", tb.tb_angle);
	fprintf(fp, "%.10f %.10f %.10f\n", tb.tb_axis[0], tb.tb_axis[1], tb.tb_axis[2]);
    	
	for (i = 0; i < 16; i++) fprintf(fp, "%.10f ", tb.tb_transform[i]);
	fprintf(fp, "\n");

	fprintf(fp, "%.10f\n", tb.tb_pan_x);
	fprintf(fp, "%.10f\n", tb.tb_pan_y);
	fprintf(fp, "%.10f\n", tb.tb_zoom);

	fclose(fp);
}


void tb_load(const char *file)
{
	int i;
    FILE *fp = fopen(file, "r");
	if (fp == NULL)
	{
		fprintf(fp, "TrackBall(): Cannot open file \"%s\" for reading.\n", file);
		return;
	}

	fscanf(fp, "%lf", &tb.tb_angle);
	fscanf(fp, "%lf %lf %lf", &tb.tb_axis[0], &tb.tb_axis[1], &tb.tb_axis[2]);
		
	for (i = 0; i < 16; i++) fscanf(fp, "%lf", &tb.tb_transform[i]);

	fscanf(fp, "%lf", &tb.tb_pan_x);
	fscanf(fp, "%lf", &tb.tb_pan_y);
	fscanf(fp, "%lf", &tb.tb_zoom);

	fclose(fp);
}

GLdouble tb_get_zoom(void) {
    return tb.tb_zoom_all;
}

void tb_get_pan(GLdouble pan[2]) {
    pan[0] = tb.tb_pan_x;
    pan[1] = tb.tb_pan_y;
}
