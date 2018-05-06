#ifndef OBJ_PARSER_H
#define OBJ_PARSER_H

#include "vec3.h"

#include <string>
#include <vector>
using namespace std;

namespace Renzoku {

#define OBJ_FILENAME_LENGTH 256
#define MATERIAL_NAME_SIZE 256
#define OBJ_LINE_SIZE 512
#define MAX_VERTEX_COUNT 4 //can only handle quads or triangles

struct obj_face
{
	int vertex_index[MAX_VERTEX_COUNT];
	int normal_index[MAX_VERTEX_COUNT];
	int texture_index[MAX_VERTEX_COUNT];
	int vertex_count;
	int material_index;
    int smooth_group;
    int object_index;
};

struct obj_sphere
{
	int pos_index;
	int up_normal_index;
	int equator_normal_index;
	int texture_index[MAX_VERTEX_COUNT];
	int material_index;
};

struct obj_plane
{
	int pos_index;
	int normal_index;
	int rotation_normal_index;
	int texture_index[MAX_VERTEX_COUNT];
	int material_index;
};

struct obj_material
{
	char name[MATERIAL_NAME_SIZE];
	char texture_filename[OBJ_FILENAME_LENGTH];
    char texture_displacement_filename[OBJ_FILENAME_LENGTH];
	float amb[3];
	float diff[3];
	float spec[3];
    float emitt[3];
	float reflect;
	float refract;
	float trans;
	float shiny[2];
	float glossy;
	float refract_index;
};

struct obj_camera
{
	int camera_pos_index;
	int camera_look_point_index;
	int camera_up_norm_index;
    float camera_focal;
    float camera_near_plane;
    float camera_far_plane;
};

struct obj_light_point
{
	int pos_index;
	int material_index;
};

struct obj_light_disc
{
	int pos_index;
	int normal_index;
	int material_index;
};

struct obj_light_quad
{
	int vertex_index[MAX_VERTEX_COUNT];
	int material_index;
};

struct obj_scene_data
{
	char scene_filename[OBJ_FILENAME_LENGTH];
	char material_filename[OBJ_FILENAME_LENGTH];
	
	vector<Vec3> vertex_list;
	vector<Vec3> vertex_normal_list;
	vector<Vec3> vertex_texture_list;
	
	vector<obj_face>   face_list;
	vector<obj_sphere> sphere_list;
	vector<obj_plane>  plane_list;
	
	vector<obj_light_point> light_point_list;
	vector<obj_light_quad> light_quad_list;
	vector<obj_light_disc> light_disc_list;
	
	vector<obj_material> material_list;
	
	vector<obj_camera> camera_list;

    vector<string> object_names;

    vector<int> smooth_groups;
};

int parse_obj_scene(obj_scene_data *data_out, const char *filename);

} // end namespace

#endif
