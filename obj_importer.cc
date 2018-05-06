#include "obj_importer.h"

#include "camera.h"
#include "triangle.h"
#include "sphere.h"
#include "quad.h"
#include "area_light.h"
#include "point_light.h"
#include "material.h"
#include "textured_material.h"
#include "bsdf.h"
#include "lambertian.h"
#include "mirror.h"
#include "glass.h"
#include "thin_transparent.h"
#include "phong.h"
#include "ward.h"
#include "surface.h"
#include "aggregate.h"
#include "aggregate_bvh.h"
#include "scene.h"
#include "log.h"

namespace Renzoku {

ObjImporter::ObjImporter() {
    
}

ObjImporter::~ObjImporter() {

}

void ObjImporter::get_shapes(Shapes &shapes) {
    shapes.assign(this->shapes->begin(), this->shapes->end());
}

void ObjImporter::get_materials(Materials &materials) {
    materials.assign(this->materials->begin(), this->materials->end());
}

void ObjImporter::get_textures(Textures &textures) {
    textures.assign(this->textures->begin(), this->textures->end());
}

void ObjImporter::get_surfaces(Surfaces &surfaces) {
    surfaces.assign(this->surfaces->begin(), this->surfaces->end());
}

void ObjImporter::get_lights(Lights &lights) {
    lights.assign(this->lights->begin(), this->lights->end());
}

void ObjImporter::get_object_names(vector<string> &object_names) {
    object_names.assign(this->object_names.begin(), this->object_names.end());
}

Scene *ObjImporter::get_scene() {
    if (lights->size() > 0 && camera) {
        Aggregate *agg = new AggregateBvh(*surfaces);
        Random *rd = new Random();
        this->scene = new Scene(name, agg, shapes, materials, textures, surfaces, lights, camera, rd);
        return scene;
    }
    return NULL;
}

static bool is_area_light(obj_material *om) {    
    if (strstr(om->name, "rz_light") != NULL) {
        //|| strstr(om->name, "emitter") != NULL || 
		//strstr(om->name, "light") != NULL)

        // light or emitter is too common name that should never be used. If a surface is misregarded as a light,
        // it will be rendered as black for indirect illumination.
        return true;

    }
    return false;
}

static void load_material(obj_material *om, string scene_path,
                          Materials *materials, Textures *textures) {

    Rgb Kd(om->diff[0], om->diff[1], om->diff[2]);
    Rgb Ks(om->spec[0], om->spec[1], om->spec[2]);
    Rgb Ke(om->emitt[0], om->emitt[1], om->emitt[2]);
    
    // HACK: don't normalize if this material represents a light
    if (is_area_light(om) == false) {
        Rgb total = Kd + Ks;
        if (total.red() > 1.0f || 
            total.green() > 1.0f || 
            total.blue() > 1.0f) {

            Log::warn() << "Glossy BRDF Kd Ks adjusted to ensure Kd + Ks <= 1." << endn;

            Kd = Kd / total.max();
            Ks = Ks / total.max();
        }
    }

    // NOTE: need to pass Ke in
    Bsdf *bsdf;
    Material *mtl;

    if (om->refract_index > 1) {
        if (om->refract_index < 1.01f) {
            // since it is very near to pass through, treat it as a pass-through material. 
            bsdf = new ThinTransparent(Kd);
        } else {
            // glass
            bsdf = new Glass(om->refract_index, Kd);
            Log::info() << "Glass loaded." << endn;
        }
    } else if (Kd.max() != (Float)0. && Ks.max() == (Float)0.) {
        // Lambertian
        bsdf = new Lambertian(Kd.red(), Kd.green(), Kd.blue());
    }         
    /*    
    else if(Kd.max() == 0. && Ks.max() != 0.) {
        mtl = new Mirror();
    } */    
    else {
        // default glossy material in OBJ will load Phong
        if (om->shiny[1] == 0.0f) {
            if (om->shiny[0] == 0.0f) {
                Log::warn() << "Phong BRDF: glossy coefficient zero. Auto set to 1." << endn;
                om->shiny[0] = 1.0f;
            }
            bsdf = new ModifiedPhong(Kd, Ks, om->shiny[0]);
        } else {
            if (om->shiny[0] == 0.0f) {
                Log::warn() << "Ward BRDF: glossy coefficient zero. Auto set to " << om->shiny[1] << "." << endn;
                om->shiny[0] = om->shiny[1];
            }
            bsdf = new Ward(Kd, Ks, 1.f / om->shiny[0], 1.f / om->shiny[1]);            
        }
    }

    // texture
    Texture *tex = NULL;
    Texture *tex_disp = NULL;
    if (om->texture_filename[0] != '\0') {
        string texture_filepath = om->texture_filename;
        for (int i = 0; i < texture_filepath.size(); ++i)
            if (texture_filepath[i] == '\\')
                texture_filepath[i] = '/';

        texture_filepath = scene_path + texture_filepath;
        tex = new Texture(texture_filepath);
    }

    if (om->texture_displacement_filename[0] != '\0') {
        string texture_filepath = om->texture_displacement_filename;
        for (int i = 0; i < texture_filepath.size(); ++i)
            if (texture_filepath[i] == '\\')
                texture_filepath[i] = '/';
        
        texture_filepath = scene_path + texture_filepath;
        tex_disp = new Texture(texture_filepath);
    }

    // TODO: to support displacement
    if (tex) {
        textures->push_back(tex);
        mtl = new TexturedMaterial(new BsdfMaterial(bsdf), tex);
    } else
        mtl = new BsdfMaterial(bsdf);
    mtl->set_name(om->name);

    materials->push_back(mtl);
}

void ObjImporter::set_home(const string home) {
    scene_path = home;
}

static void export_mitsuba(Materials &materials) {
    for (int i = 0; i < materials.size(); ++i) {
        Material *m = materials[i];
        
        switch (m->get_bsdf()->get_bsdf_type()) {
            case Bsdf::LAMBERTIAN:
            {
                Lambertian *lambert = (Lambertian *)m->get_bsdf();
                Rgb diffuse = lambert->get_diffuse_component();            
                cout << "<bsdf id=\"" << m->get_name() << "_material\" type=\"diffuse\">" << endl;
                cout << "    <rgb name=\"reflectance\" value=\"" << diffuse.red() << " " << diffuse.green() << " " << diffuse.blue() << "\"/>" << endl;            
                cout << "</bsdf>" << endl;
                break;
            }
            case Bsdf::PHONG:
            {
                ModifiedPhong *phong = (ModifiedPhong *)m->get_bsdf();
                Rgb diffuse = phong->get_diffuse_component();
                Rgb specular = phong->get_specular_component();
                cout << "<bsdf id=\"" << m->get_name() << "_material\" type=\"phong\">" << endl;
                cout << "    <rgb name=\"diffuseReflectance\" value=\"" << diffuse.red() << " " << diffuse.green() << " " << diffuse.blue() << "\"/>" << endl;
                cout << "    <rgb name=\"specularReflectance\" value=\"" << specular.red() << " " << specular.green() << " " << specular.blue() << "\"/>" << endl;
                cout << "    <float name=\"exponent\" value=\"" << phong->get_gloss_exponential() << "\"/>" << endl;
                cout << "</bsdf>" << endl;
                break;
            }

            default:
                Log::info() << "Mitsuba equivalent material not yet supported." << endn;
                break;
        }
    }
}

void ObjImporter::load(const string file) {
    string name = file;
    unsigned int found = name.find_last_of('.');
    name = name.substr(0, found);

    int loaded = parse_obj_scene(&obj_data, file.c_str());
	if(! loaded)
	{
        Log::critical("OBJ parse error.", ExitCode::EXIT_CODE_SCENE_FILE_ERROR);
    }
    
    Materials   *materials  = new Materials;
    Shapes      *shapes     = new Shapes;
    Surfaces    *surfaces   = new Surfaces;
    Lights      *lights     = new Lights;
    Textures    *textures   = new Textures;
    vector<bool> is_light(obj_data.material_list.size());
    Camera *camera          = NULL;

    // customization to add camera and light to OBJ file
    {
        //---Import the camera        
        if (obj_data.camera_list.size() > 0) {
			obj_camera *c = &obj_data.camera_list[0];
            camera = new Camera();
            Sensor *sensor = new Sensor(Size2(512, 512), Size2(0.025, 0.025));
     
			Vec3 pos = obj_data.vertex_list[c->camera_pos_index];
          	Vec3 lookat = obj_data.vertex_list[c->camera_look_point_index];
         	Vec3 upn = obj_data.vertex_normal_list[c->camera_up_norm_index];
            Float focal = c->camera_focal;    
            camera->set_perspective(pos, lookat, upn, focal, sensor);
            //camera->set(pos, lookat, upn, focal, Size2(0.0375, 0.025), Size2(768, 512));
            if (c->camera_near_plane > 0) {
                camera->set_near_plane(c->camera_near_plane);
            }
            if (c->camera_far_plane > 0) {
                camera->set_far_plane(c->camera_far_plane);
            }
        }

        //---Import the light sources
        //@1 Point lights
        for(int i=0; i<obj_data.light_point_list.size(); i++) {
            //Point light position
            obj_light_point *lp = &obj_data.light_point_list[i];
            Vec3 p = obj_data.vertex_list[lp->pos_index];
            Light *pl = new PointLight(p, DefaultRgb::white);
            lights->push_back(pl);
        }
        //@2 Area Lights
        for(int i=0; i<obj_data.light_quad_list.size(); i++) {
            obj_light_quad *lq = &obj_data.light_quad_list[i];
           
            Vec3 p0 = obj_data.vertex_list[lq->vertex_index[0]];
            Vec3 p1 = obj_data.vertex_list[lq->vertex_index[1]];
            Vec3 p2 = obj_data.vertex_list[lq->vertex_index[2]];
            Vec3 p3 = obj_data.vertex_list[lq->vertex_index[3]];

            obj_material *om = &obj_data.material_list[lq->material_index];        
            Surface al(new Quad(p0, p1, p2, p3), Rgb(om->diff[0], om->diff[1], om->diff[2]));

            lights->push_back(al.get_area_light());
            surfaces->append(al);
        }
    }
    
    //---Import the materials    
    for(int i = 0; i < obj_data.material_list.size(); ++i) {
    	obj_material *om = &obj_data.material_list[i];

        // HACK: use a special material name to indicate area light
        is_light[i] = is_area_light(om);

        load_material(om, scene_path, materials, textures);
    }
    
    //export_mitsuba(*materials);

    //---Import object names    
    object_names.assign(obj_data.object_names.begin(),
                        obj_data.object_names.end());
    

    //---Import the surfaces
    //@1 Push the faces into the queue
    for(int i=0; i<obj_data.face_list.size(); i++) {
        obj_face *f = &obj_data.face_list[i];                
        
        Shape *s = NULL;
        if (f->vertex_count == 3) {            
            Vec3 p0 = obj_data.vertex_list[f->vertex_index[0]];
            Vec3 p1 = obj_data.vertex_list[f->vertex_index[1]];
            Vec3 p2 = obj_data.vertex_list[f->vertex_index[2]];
        
            // texture coordinates
            Vec2 uv0, uv1, uv2;
            if (f->texture_index[0] < 0) {
                //Log::info() << "Invalid texture index: " << f->texture_index[0] << endn;
            } else {
                Vec3 v = obj_data.vertex_texture_list[f->texture_index[0]];
                uv0 = Vec2(v.e[0], v.e[1]);
            }

            if (f->texture_index[1] < 0) {
                //Log::info() << "Invalid texture index: " << f->texture_index[1] << endn;
            } else {            
                Vec3 v = obj_data.vertex_texture_list[f->texture_index[1]];
                uv1 = Vec2(v.e[0], v.e[1]);
            }

            if (f->texture_index[2] < 0) {
                //Log::info() << "Invalid texture index: " << f->texture_index[2] << endn;
            } else {            
                Vec3 v = obj_data.vertex_texture_list[f->texture_index[2]];
                uv2 = Vec2(v.e[0], v.e[1]);
            }
            
            /*
            if (uv0.x() < 0 || uv0.y() < 0 || uv1.x() < 0 || uv1.y() < 0 ||
                uv2.x() < 0 || uv2.y() < 0) {
                Log::warn() << "Negative UV." << endn;
            }*/

            s = new Triangle(p0, p1, p2, uv0, uv1, uv2);
                        
            if (f->normal_index[0] >= 0) {
                s->set_vertex_normal(0, obj_data.vertex_normal_list[f->normal_index[0]]);
            }
            if (f->normal_index[1] >= 0) {
                s->set_vertex_normal(1, obj_data.vertex_normal_list[f->normal_index[1]]);
            }
            if (f->normal_index[2] >= 0) {
                s->set_vertex_normal(2, obj_data.vertex_normal_list[f->normal_index[2]]);
            }
            
        } else if (f->vertex_count == 4) {
            Vec3 p0(obj_data.vertex_list[f->vertex_index[0]]);
            Vec3 p1(obj_data.vertex_list[f->vertex_index[1]]);
            Vec3 p2(obj_data.vertex_list[f->vertex_index[2]]);
            Vec3 p3(obj_data.vertex_list[f->vertex_index[3]]);

            // texture coordinates
            Vec2 uv0, uv1, uv2, uv3;
            if (f->texture_index[0] < 0) {
                //Log::info() << "Invalid texture index: " << f->texture_index[0] << endn;
            } else {
                Vec3 v = obj_data.vertex_texture_list[f->texture_index[0]];
                uv0 = Vec2(v.e[0], v.e[1]);
            }

            if (f->texture_index[1] < 0) {
                //Log::info() << "Invalid texture index: " << f->texture_index[1] << endn;
            } else {            
                Vec3 v = obj_data.vertex_texture_list[f->texture_index[1]];
                uv1 = Vec2(v.e[0], v.e[1]);
            }

            if (f->texture_index[2] < 0) {
                //Log::info() << "Invalid texture index: " << f->texture_index[2] << endn;
            } else {            
                Vec3 v = obj_data.vertex_texture_list[f->texture_index[2]];
                uv2 = Vec2(v.e[0], v.e[1]);
            }
            
            if (f->texture_index[3] < 0) {
                //Log::info() << "Invalid texture index: " << f->texture_index[3] << endn;
            } else {            
                Vec3 v = obj_data.vertex_texture_list[f->texture_index[3]];
                uv3 = Vec2(v.e[0], v.e[1]);
            }
            /*
            if (uv0.x() < 0 || uv0.y() < 0 || uv1.x() < 0 || uv1.y() < 0 ||
                uv2.x() < 0 || uv2.y() < 0 || uv3.x() < 0 || uv3.y() < 0) {
                Log::warn() << "Negative UV." << endn;
            }*/

            s = new Quad(p0, p1, p2, p3, uv0, uv1, uv2, uv3);

            if (f->normal_index[0] >= 0) {
                s->set_vertex_normal(0, obj_data.vertex_normal_list[f->normal_index[0]]);
            }
            if (f->normal_index[1] >= 0) {
                s->set_vertex_normal(1, obj_data.vertex_normal_list[f->normal_index[1]]);
            }
            if (f->normal_index[2] >= 0) {
                s->set_vertex_normal(2, obj_data.vertex_normal_list[f->normal_index[2]]);
            }
            if (f->normal_index[3] >= 0) {
                s->set_vertex_normal(3, obj_data.vertex_normal_list[f->normal_index[3]]);
            }
            
        } else {
            // FIXME: use country_kitchen scene and fix this.
            //Log::info() << "Polygon data: " << obj_data.vertex_list[f->vertex_index[0]] << endn;
            printf("Polygon with vertex count %d not supported. Surface skipped.\n", f->vertex_count);            
            continue;
        }

        // OBJ format: 0 means no group
        // in Renzoku, group number starts from zero.
        // only let shape participate in smoothing when no vertex normal is explicitly given.
        // assume if one vertex normal is not given, other vertices in the shape are not either.
        if (f->smooth_group > 0 && f->normal_index[0] < 0)
            s->set_smooth_group(f->smooth_group);

        s->set_object_index(f->object_index);

        shapes->push_back(s);
        
        if (f->material_index >= 0) {
            if (is_light[f->material_index]) {
                Material *m = materials->at(f->material_index);
                IDiffuse *bsdf = dynamic_cast<IDiffuse *>(m->get_bsdf());
                if (bsdf) {
                    // TODO: extend to textured area light
                    Rgb kd = bsdf->get_diffuse_component();
                    Surface al(s, kd);
                    lights->push_back(al.get_area_light());
                    surfaces->append(al);
                } else {
                    Log::warn() << "Area light specification invalid." << endn;
                }
                // TODO: the dummy material m should be deleted. 
            } else {
                Material *m = materials->at(f->material_index);        
                Surface sf(s, m);
                surfaces->append(sf);
            }
        }
    }

    //@2 Push the spheres into the queue
    for(int i=0; i<obj_data.sphere_list.size(); i++) {
        obj_sphere *s = &obj_data.sphere_list[i];
        Vec3 c(obj_data.vertex_list[s->pos_index]);//sphere center
        Vec3 upn(obj_data.vertex_list[s->up_normal_index]);
        //The length of the up normal is the radius
        Shape *temp = new Sphere(c, upn.length());

        //~~~Connect the geometry with its material property
        Material *m = materials->at(s->material_index);        

        Surface sf(temp, m);
        surfaces->append(sf);
    }

    // smooth groups vertex normal
    if (obj_data.smooth_groups.size() > 0) {
        Log::info() << "Computing vertex normals for smooth groups..." << endn;
                
        vector<Vec3> vertex_list(4);
        vector<Vec3> vertex_list2(4);
        for (int g = 0; g < obj_data.smooth_groups.size(); ++g) {
            int group = obj_data.smooth_groups[g];
            Shapes shape_group;

            int count = 0;
            for (int i = 0; i < shapes->size(); ++i) {
                Shape *s = (*shapes)[i];
                if (s->get_smooth_group() == group) count++;
            }

            shape_group.resize(count);
            count = 0;
            for (int i = 0; i < shapes->size(); ++i) {
                Shape *s = (*shapes)[i];
                if (s->get_smooth_group() == group) 
                    shape_group[count++] = s;
            }

            // find normal per vertex
            for (int i = 0; i < shape_group.size(); ++i) {
                Shape *s = shape_group[i];
            
                vertex_list.clear();
                s->get_vertex_positions(vertex_list);
                for (int j = 0; j < vertex_list.size(); ++j) {
                    Vec3 n;

                    for (int k = 0; k < shape_group.size(); ++k) {
                        Shape *p = shape_group[k];
                
                        vertex_list2.clear();
                        p->get_vertex_positions(vertex_list2);

                        for (int h = 0; h < vertex_list2.size(); ++h) {
                            if (vertex_list[j] == vertex_list2[h]) {
                                n += p->normal(vertex_list[j]);
                            }
                        }
                    }

                    Vec3 normalized_n = unit_vector(n);
                    if (normalized_n.is_nan())
                        s->set_vertex_normal(j, Vec3(0.0f));
                    else
                        s->set_vertex_normal(j, normalized_n);
                }
            }       
        }

        /*
        for (int i = 0; i < shapes->size(); ++i) {
            Shape *s = (*shapes)[i];
            if (s->get_smooth_group() < 0) continue;

            vector<Vec3> vertex_list;
            s->get_vertex_positions(vertex_list);
            for (int j = 0; j < vertex_list.size(); ++j) {
                Vec3 n;

                for (int k = 0; k < shapes->size(); ++k) {
                    Shape *p = (*shapes)[k];
                    if (p->get_smooth_group() == s->get_smooth_group()) {
                        vector<Vec3> vertex_list2;
                        p->get_vertex_positions(vertex_list2);

                        for (int h = 0; h < vertex_list2.size(); ++h) {
                            if (vertex_list[j] == vertex_list2[h]) {
                                n += p->normal(vertex_list[j]);
                            }
                        }
                    }
                }

                s->set_vertex_normal(j, unit_vector(n));
            }        
        }
        */
    }

    // assign results
    this->name = name;
    this->shapes = shapes;
    this->materials = materials;
    this->textures = textures;
    this->surfaces = surfaces;
    this->lights = lights;
    this->camera = camera;
}

} // end namespace