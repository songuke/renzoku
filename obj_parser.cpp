#include "obj_parser.h"
#include "log.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

namespace Renzoku {

#define WHITESPACE " \t\n\r"

static char strequal(const char *s1, const char *s2)
{
	if (strcmp(s1, s2) == 0)
		return 1;
	return 0;
}

static char contains(const char *haystack, const char *needle)
{
	if (strstr(haystack, needle) == NULL)
		return 0;
	return 1;
}

static int obj_convert_to_list_index(int current_max, int index)
{
	if(index == 0)  //no index
		return -1;
		
	if(index < 0)  //relative to current list position
		return current_max + index;
		
	return index - 1;  //normal counting index
}

static void obj_convert_to_list_index_v(int current_max, int *indices)
{
	for(int i=0; i<MAX_VERTEX_COUNT; i++)
		indices[i] = obj_convert_to_list_index(current_max, indices[i]);
}

static void obj_set_material_defaults(obj_material *mtl)
{
	mtl->amb[0] = 0.0;
	mtl->amb[1] = 0.0;
	mtl->amb[2] = 0.0;
	mtl->diff[0] = 0.0;
	mtl->diff[1] = 0.0;
	mtl->diff[2] = 0.0;
	mtl->spec[0] = 0.0;
	mtl->spec[1] = 0.0;
	mtl->spec[2] = 0.0;
    mtl->emitt[0] = 0.0;
    mtl->emitt[1] = 0.0;
    mtl->emitt[2] = 0.0;
	mtl->reflect = 0.0;
	mtl->trans = 1;
	mtl->glossy = 1;
	mtl->shiny[0] = 0.0;
    mtl->shiny[1] = 0.0;
	mtl->refract_index = 1;
	mtl->texture_filename[0] = '\0';
    mtl->texture_displacement_filename[0] = '\0';
}

static int obj_parse_vertex_index(int *vertex_index, int *texture_index, int *normal_index)
{
	char *temp_str;
	char *token;
	int vertex_count = 0;

	
	while( (token = strtok(NULL, WHITESPACE)) != NULL)
	{
        if (vertex_count >= MAX_VERTEX_COUNT) {
            printf("Polygon with vertices more than %d are not supported. Vertex index from %d is ignored.\n", MAX_VERTEX_COUNT, MAX_VERTEX_COUNT);
            // TODO: to support more than 4 vertices surfaces.
            // must continue otherwise stack might corrupt.
            continue;
        }

		if(texture_index != NULL)
			texture_index[vertex_count] = 0;
		if(normal_index != NULL)
		    normal_index[vertex_count] = 0;

		vertex_index[vertex_count] = atoi( token );
		
		if(contains(token, "//"))  //normal only
		{
			temp_str = strchr(token, '/');
			temp_str++;
			normal_index[vertex_count] = atoi( ++temp_str );
		}
		else if(contains(token, "/"))
		{
			temp_str = strchr(token, '/');
			texture_index[vertex_count] = atoi( ++temp_str );

			if(contains(temp_str, "/"))
			{
				temp_str = strchr(temp_str, '/');
				normal_index[vertex_count] = atoi( ++temp_str );
			}
		} else {
            // only vertex index, no other info.
        }
		
		vertex_count++;
	}

	return vertex_count;
}

static void obj_parse_face(obj_scene_data *scene, obj_face* face)
{
	int vertex_count;
	
	vertex_count = obj_parse_vertex_index(face->vertex_index, face->texture_index, face->normal_index);
	obj_convert_to_list_index_v(scene->vertex_list.size(), face->vertex_index);
	obj_convert_to_list_index_v(scene->vertex_texture_list.size(), face->texture_index);
	obj_convert_to_list_index_v(scene->vertex_normal_list.size(), face->normal_index);
	face->vertex_count = vertex_count;
}

static void obj_parse_sphere(obj_scene_data *scene, obj_sphere *sphere)
{
	int temp_indices[MAX_VERTEX_COUNT];

	obj_parse_vertex_index(temp_indices, sphere->texture_index, NULL);
	obj_convert_to_list_index_v(scene->vertex_texture_list.size(), sphere->texture_index);
	sphere->pos_index = obj_convert_to_list_index(scene->vertex_list.size(), temp_indices[0]);
	sphere->up_normal_index = obj_convert_to_list_index(scene->vertex_normal_list.size(), temp_indices[1]);
	sphere->equator_normal_index = obj_convert_to_list_index(scene->vertex_normal_list.size(), temp_indices[2]);
}

static void obj_parse_plane(obj_scene_data *scene, obj_plane *plane)
{
	int temp_indices[MAX_VERTEX_COUNT];

	obj_parse_vertex_index(temp_indices, plane->texture_index, NULL);
	obj_convert_to_list_index_v(scene->vertex_texture_list.size(), plane->texture_index);
	plane->pos_index = obj_convert_to_list_index(scene->vertex_list.size(), temp_indices[0]);
	plane->normal_index = obj_convert_to_list_index(scene->vertex_normal_list.size(), temp_indices[1]);
	plane->rotation_normal_index = obj_convert_to_list_index(scene->vertex_normal_list.size(), temp_indices[2]);
}

static void obj_parse_light_point(obj_scene_data *scene, obj_light_point *o)
{
	o->pos_index = obj_convert_to_list_index(scene->vertex_list.size(), atoi( strtok(NULL, WHITESPACE)) );
}

static void obj_parse_light_quad(obj_scene_data *scene, obj_light_quad *o)
{
	obj_parse_vertex_index(o->vertex_index, NULL, NULL);
	obj_convert_to_list_index_v(scene->vertex_list.size(), o->vertex_index);
}

static void obj_parse_light_disc(obj_scene_data *scene, obj_light_disc *o)
{
	int temp_indices[MAX_VERTEX_COUNT];

	obj_parse_vertex_index(temp_indices, NULL, NULL);
	o->pos_index = obj_convert_to_list_index(scene->vertex_list.size(), temp_indices[0]);
	o->normal_index = obj_convert_to_list_index(scene->vertex_normal_list.size(), temp_indices[1]);
}

Vec3 obj_parse_vector()
{
	Vec3 v;
	v.e[0] = atof( strtok(NULL, WHITESPACE));
	v.e[1] = atof( strtok(NULL, WHITESPACE));
	v.e[2] = atof( strtok(NULL, WHITESPACE));
	return v;
}

Vec3 obj_parse_vector2or3()
{
	Vec3 v;
	v.e[0] = atof( strtok(NULL, WHITESPACE));
	v.e[1] = atof( strtok(NULL, WHITESPACE));
    
    // some vector (such as vertex texture coordinates, are defined with 2 or 3 elements)
    char *token;
    if ( (token = strtok(NULL, WHITESPACE)) != NULL ) {
	    v.e[2] = atof( token );
    } else {
        v.e[2] = 0.0f;
    }
	return v;
}

static void obj_parse_camera(obj_scene_data *scene, obj_camera *camera)
{
	int indices[3];
	int count = 0;

    char *token;	
	while( (token = strtok(NULL, WHITESPACE)) != NULL )
	{
	    indices[count] = atoi( token );
        count++;

        if (count == 3) break;
	}
    
	camera->camera_pos_index = obj_convert_to_list_index(scene->vertex_list.size(), indices[0]);
	camera->camera_look_point_index = obj_convert_to_list_index(scene->vertex_list.size(), indices[1]);
	camera->camera_up_norm_index = obj_convert_to_list_index(scene->vertex_normal_list.size(), indices[2]);
    camera->camera_focal = atof( strtok(NULL, WHITESPACE) );

    if ((token = strtok(NULL, WHITESPACE)) != NULL) {
        camera->camera_near_plane = atof( token );
    }
    if ((token = strtok(NULL, WHITESPACE)) != NULL) {
        camera->camera_far_plane = atof( token );
    }

    // TODO: support film size, image size here
}

int obj_parse_mtl_file(const char *filename, vector<obj_material> &material_list)
{
	int line_number = 0;
	char *current_token;
	char current_line[OBJ_LINE_SIZE];
	char material_open = 0;
	obj_material *current_mtl = NULL;
	FILE *mtl_file_stream;

	// open scene
	mtl_file_stream = fopen(filename, "r");
	if (mtl_file_stream == 0)
	{
		fprintf(stderr, "Error reading file: %s\n", filename);
		return 0;
	}

	while (fgets(current_line, OBJ_LINE_SIZE, mtl_file_stream))
	{
		current_token = strtok(current_line, WHITESPACE);
		line_number++;

		//skip comments
		if (current_token == NULL || strequal(current_token, "//") || strequal(current_token, "#"))
			continue;


		//start material
		else if (strequal(current_token, "newmtl"))
		{
			material_open = 1;

			if (current_mtl)
				material_list.push_back(*current_mtl);
			else
				current_mtl = new obj_material;

			// make new material
			obj_set_material_defaults(current_mtl);

			// get the name
			strncpy(current_mtl->name, strtok(NULL, WHITESPACE), MATERIAL_NAME_SIZE);
		}

		//ambient
		else if (strequal(current_token, "Ka") && material_open)
		{
			current_mtl->amb[0] = atof(strtok(NULL, " \t"));
			current_mtl->amb[1] = atof(strtok(NULL, " \t"));
			current_mtl->amb[2] = atof(strtok(NULL, " \t"));
		}

		//diff
		else if (strequal(current_token, "Kd") && material_open)
		{
			current_mtl->diff[0] = atof(strtok(NULL, " \t"));
			current_mtl->diff[1] = atof(strtok(NULL, " \t"));
			current_mtl->diff[2] = atof(strtok(NULL, " \t"));
		}

		//specular
		else if (strequal(current_token, "Ks") && material_open)
		{
			current_mtl->spec[0] = atof(strtok(NULL, " \t"));
			current_mtl->spec[1] = atof(strtok(NULL, " \t"));
			current_mtl->spec[2] = atof(strtok(NULL, " \t"));
		}

		// emission
		else if (strequal(current_token, "Ke") && material_open)
		{
			current_mtl->emitt[0] = atof(strtok(NULL, " \t"));
			current_mtl->emitt[1] = atof(strtok(NULL, " \t"));
			current_mtl->emitt[2] = atof(strtok(NULL, " \t"));
		}

		//shiny
		else if (strequal(current_token, "Ns") && material_open)
		{

			current_mtl->shiny[0] = atof(strtok(NULL, " \t"));
			if (current_mtl->shiny[0] == 0.0f) {
				Log::warn() << "Ns parameter of material " << current_mtl->name << " is zero." << endn;
			}

			// EXTEND: support anisotropic glossy for Ward model
			char *token;
			if ((token = strtok(NULL, " \t")) != NULL) {
				current_mtl->shiny[1] = atof(token);
				if (current_mtl->shiny[1] == 0.0f) {
					Log::warn() << "Ns parameter of material " << current_mtl->name << " is zero." << endn;
				}
			}
			else {
				// let shiny[1] not assigned to differentiate between Phong and Ward
			}
		}
		//transparent
		else if (strequal(current_token, "d") && material_open)
		{
			current_mtl->trans = atof(strtok(NULL, " \t"));
		}
		//reflection
		else if (strequal(current_token, "r") && material_open)
		{
			current_mtl->reflect = atof(strtok(NULL, " \t"));
		}
		//glossy
		else if (strequal(current_token, "sharpness") && material_open)
		{
			current_mtl->glossy = atof(strtok(NULL, " \t"));
		}
		//refract index
		else if (strequal(current_token, "Ni") && material_open)
		{
			current_mtl->refract_index = atof(strtok(NULL, " \t"));
		}
		// illumination type
		else if (strequal(current_token, "illum") && material_open)
		{
		}
		// texture map
		else if (strequal(current_token, "map_Kd") && material_open)
		{
			strncpy(current_mtl->texture_filename, strtok(NULL, WHITESPACE), OBJ_FILENAME_LENGTH);
		}
		else if (strequal(current_token, "map_Disp") && material_open)
		{
			strncpy(current_mtl->texture_displacement_filename, strtok(NULL, WHITESPACE), OBJ_FILENAME_LENGTH);
		}
		else
		{
			fprintf(stderr, "Unknown command '%s' in material file %s at line %i:\n\t%s\n",
				current_token, filename, line_number, current_line);
			//return 0;
		}
	}

	if (current_mtl) {
		material_list.push_back(*current_mtl);
		delete current_mtl;
	}

	fclose(mtl_file_stream);

	return 1;

}

int obj_parse_obj_file(obj_scene_data *growable_data, const char *filename)
{
	FILE* obj_file_stream;	
	char *current_token = NULL;
	char current_line[OBJ_LINE_SIZE];
	int line_number = 0;

	obj_file_stream = fopen( filename, "r");
	if(obj_file_stream == 0)	// failed to open file
	{
		fprintf(stderr, "Error reading file: %s\n", filename);
		return 0;
	}

    int current_material     = -1;		 // default no material
    
    int current_smooth_group = 0;        // zero no vertex smoothing        
    
    int current_object_index = -1;       // default no object specified

    // pre-allocate some memory to avoid repeat expanding for large scenes
    int prealloc = 1024 * 1024;
    growable_data->vertex_list.reserve(prealloc);
    growable_data->vertex_normal_list.reserve(prealloc);
    growable_data->vertex_texture_list.reserve(prealloc);
    growable_data->face_list.reserve(prealloc);

	// parser loop
	while( fgets(current_line, OBJ_LINE_SIZE, obj_file_stream) )
	{
		current_token = strtok( current_line, WHITESPACE);
		line_number++;
		
		// skip comments
		if( current_token == NULL || current_token[0] == '#')
			continue;

		// parse objects
		else if( strequal(current_token, "v") )  // process vertex
		{
			growable_data->vertex_list.push_back( obj_parse_vector() );
		}
		
		else if( strequal(current_token, "vn") )  // process vertex normal
		{
			growable_data->vertex_normal_list.push_back( obj_parse_vector() );
		}
		
		else if( strequal(current_token, "vt") )  // process vertex texture
		{
			growable_data->vertex_texture_list.push_back( obj_parse_vector2or3() );
		}
		
		else if( strequal(current_token, "f") )  // process face
		{
			obj_face face;
			obj_parse_face(growable_data, &face);
			face.material_index = current_material;
            face.object_index = current_object_index;
            face.smooth_group = current_smooth_group;
			growable_data->face_list.push_back(face);
		}
		
		else if( strequal(current_token, "sp") )  // process sphere
		{
			obj_sphere sphere;
			obj_parse_sphere(growable_data, &sphere);
			sphere.material_index = current_material;
			growable_data->sphere_list.push_back(sphere); 
		}
		
		else if( strequal(current_token, "pl") )  // process plane
		{
			obj_plane plane;
			obj_parse_plane(growable_data, &plane);
			plane.material_index = current_material;
			growable_data->plane_list.push_back(plane);
		}
		
		else if( strequal(current_token, "p") )  // process point
		{
			//make a small sphere to represent the point?
		}
		
		else if( strequal(current_token, "lp") )  // light point source
		{
			obj_light_point o;
			obj_parse_light_point(growable_data, &o);
			o.material_index = current_material;
			growable_data->light_point_list.push_back(o);
		}
		
		else if( strequal(current_token, "ld") )  // process light disc
		{
			obj_light_disc o;
			obj_parse_light_disc(growable_data, &o);
			o.material_index = current_material;
			growable_data->light_disc_list.push_back(o);
		}
		
		else if( strequal(current_token, "lq") )  // process light quad
		{
			obj_light_quad o;
			obj_parse_light_quad(growable_data, &o);
			o.material_index = current_material;
			growable_data->light_quad_list.push_back(o);
		}
		
		else if( strequal(current_token, "c") )  // camera
		{
			obj_camera camera;
			obj_parse_camera(growable_data, &camera);
			growable_data->camera_list.push_back(camera);
		}
		
		else if( strequal(current_token, "usemtl") )  // usemtl
		{ 
			bool found = false;
			string name = strtok(NULL, WHITESPACE);
			for (int i = 0; i < growable_data->material_list.size(); ++i) {
				if (name == growable_data->material_list[i].name) {
					current_material = i;
					found = true;
					break;
				}
			}

			if (! found) {
				current_material = -1;
				Log::error() << "Material " << name << " not found." << endn;
			}
		}
		
		else if( strequal(current_token, "mtllib") )  // mtllib
		{
			strncpy(growable_data->material_filename, strtok(NULL, WHITESPACE), OBJ_FILENAME_LENGTH);
            
            char mtl_file[256]; 
            for (int i = 0; i < 256; ++i) mtl_file[i] = 0;

            char *slash = (char *)filename + strlen(filename);
            do {
                slash--;
            } while (slash >= filename && *slash != '/' && *slash != '\\');
            if (slash < filename) 
                obj_parse_mtl_file(growable_data->material_filename, growable_data->material_list);
			else {
                strncpy(mtl_file, filename, slash - filename + 1);
                strcat(mtl_file, growable_data->material_filename);
                cout << "MTL file: " << mtl_file << endl;
                obj_parse_mtl_file(mtl_file, growable_data->material_list);
            }
            continue;
		}
		
		else if( strequal(current_token, "o") )  // object name
		{ 
            char *obj_name = strtok(NULL, WHITESPACE);
            growable_data->object_names.push_back(obj_name);
            current_object_index++;
        }
		else if( strequal(current_token, "s") )  // smoothing
		{ 
            char *group = strtok(NULL, WHITESPACE);
            if (strequal(group, "off")) {
                current_smooth_group = 0;
            } else {
                current_smooth_group = atof(group);

                // keep track of all smooth groups
                bool found = false;
                for (int i = 0; i < growable_data->smooth_groups.size(); ++i) {
                    if (growable_data->smooth_groups[i] == current_smooth_group) {
                        found = true;
                        break;
                    }
                }
                if (! found)
                    growable_data->smooth_groups.push_back(current_smooth_group);
            }
        }
		else if( strequal(current_token, "g") )  // group
		{ }		

		else
		{
			printf("Unknown command '%s' in scene code at line %i: \"%s\".\n",
					current_token, line_number, current_line);
		}
	}

	fclose(obj_file_stream);
	
	return 1;
}

int parse_obj_scene(obj_scene_data *data_out, const char *filename)
{
	if( obj_parse_obj_file(data_out, filename) == 0 )    // failed to parse the scene file
		return 0;
	
	return 1;
}

} // end namespace
