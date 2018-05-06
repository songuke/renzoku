#include "json_importer.h"
#include "obj_importer.h"

#include "rgb.h"
#include "lambertian.h"
#include "ward.h"
#include "phong.h"
#include "glass.h"
#include "mirror.h"
#include "thin_transparent.h"

#include "light.h"
#include "area_light.h"
#include "point_light.h"
#include "spot_light.h"
#include "dir_light.h"
#include "surface.h"

#include "camera.h"
#include "camera_view.h"

#include "aggregate_bvh.h"

#include "file.h"

#include "log.h"
#include "stats.h"

#include "monte_carlo_integrator.h";
#include "path_tracing.h"
#include "vpl.h"
//#include "vpl_gpu.h"
//#include "vsl_gpu.h"

namespace Renzoku {

string JsonImporter::keywords[] = {
    "geometry",
    "materials",
    "light_specs",
    "cameras",
    "scene",
    "world",
    "surfaces",
    "lights",
    "integrator",
    "camera",
    "type"
};

JsonImporter::JsonImporter() 
    : scene_options(NULL)
{    
}

bool JsonImporter::process_object(int stack_level, string last_key, JsonObject *obj) {
    if (stack_level == 1) {    // only takes care of special keywords in first level
        if (last_key == keywords[KW_GEOMETRY]) {
            Log::info() << "Load geometry: " << endn;
            load_geometry(obj);
            return true;
        } else
        if (last_key == keywords[KW_MATERIALS]) {
            Log::info() << "Load materials: " << endn;
            load_materials(obj);
            return true;
        } else
        if (last_key == keywords[KW_LIGHT_SPECS]) {
            Log::info() << "Load light specs: " << endn;
            load_light_specs(obj);
            return true;
        } else
        if (last_key == keywords[KW_CAMERAS]) {
            Log::info() << "Load cameras: " << endn;
            load_cameras(obj);
            return true;
        } else
        if (last_key == keywords[KW_SCENE]) {
            Log::info() << "Load scene: " << endn;
            load_scene(obj);
            return true;
        } else
        return false;
    }
    return false;
}

void JsonImporter::process_root(JsonObject *root) {
    if (root->is_empty()) {
        Log::info() << "Nothing left for further processing." << endn;
    }
}

bool JsonImporter::process_array(string last_key, JsonArray *array) {
    return false;
}

void JsonImporter::get_shape(string name, Shapes &matched) {
    matched.clear();

    for (int i = 0; i < shapes->size(); ++i) {
        string shape_name = (*shapes)[i]->get_name();
        if (shape_name.compare(name) == 0 || 
            shape_name.find(name + ".") == 0) {
        // OBJ objects might be concatenated to the geometry name defined in JSON. 
        // so when match, for example, ball, it should match with string that begins with ball. as well.

            matched.push_back((*shapes)[i]);
        }
    }
}

Material *JsonImporter::get_material(string name) {
    for (int i = 0; i < materials->size(); ++i)
        if ((*materials)[i]->get_name().compare(name) == 0)
            return (*materials)[i];

    return NULL;
}

LightSpec *JsonImporter::get_light_spec(string name) {
    for (int i = 0; i < light_specs->size(); ++i)
        if ((*light_specs)[i]->get_name().compare(name) == 0)
            return (*light_specs)[i];

    return NULL;
}

Camera *JsonImporter::get_camera(string name) {
    for (int i = 0; i < cameras->size(); ++i)
        if ((*cameras)[i]->get_name().compare(name) == 0)
            return (*cameras)[i];

    return NULL;
}

static void load_shapes_from_obj(string json_key, ObjImporter &obj_importer, Shapes &shapes) {
    obj_importer.get_shapes(shapes);
    vector<string> obj_shape_names;    // part name defined in OBJ file
    obj_importer.get_object_names(obj_shape_names);

    string name = json_key;

    // TODO: to support object group to save memory here
    for (int j = 0; j < shapes.size(); ++j) {
        int k = shapes[j]->get_object_index();
        if (k >= 0) {
            shapes[j]->set_name(name + "." + obj_shape_names[k]);
        } else {
            shapes[j]->set_name(name);
        }
    }
}

void JsonImporter::load_geometry(JsonObject *obj) {

    for (int i = 0; i < obj->keys.size(); ++i) {

        JsonValue *v = obj->values[i];
        string name = obj->keys[i];

        switch (v->type) {
        case JsonValue::JSON_VALUE_PRIMITIVE:
            if (v->data.primitive->is_string()) {
                string file = v->data.primitive->to_string();
                // only OBJ for now
                ObjImporter obj_importer;
                LogStream &out = Log::info();
                out << "Loading " << setw(30) << left << file << " ... " << endn;
                Stats stats;
                stats.tic();
                string obj_file = scene_path + "mesh/" + file;
                obj_importer.set_home(scene_path);
                obj_importer.load(obj_file);

                Shapes obj_shapes;
                load_shapes_from_obj(name, obj_importer, obj_shapes);
                shapes->insert(shapes->end(), obj_shapes.begin(), obj_shapes.end());
                stats.toc();
                out << setw(10) << " ... " << stats.elapsed() << " ms" << endn;
				
            } else {
                Log::error() << "Invalid geometry data." << endn;
            }
            break;

        case JsonValue::JSON_VALUE_ARRAY:
            Log::error() << "Unsupported geometry data." << endn;
            break;

        case JsonValue::JSON_VALUE_OBJECT:
            JsonObject* obj = v->data.object;            
            string type = obj->get_value_strict("type")->to_string();
            if (type == "sphere") {

                Vec3 center = obj->get_value_strict("center")->to_array()->to_vec3();
                Float rad = obj->get_value_strict("radius")->to_numeric();
                Sphere *s = new Sphere(center, rad);
                s->set_name(name);
                shapes->push_back(s);

            } else {

                Log::error() << "Unsupported geometry type: " << type << endn;
            
            }
            break;
        }
    }
}

void JsonImporter::load_materials(JsonObject *obj) {
    // TODO: to support texture loading with materials

    for (int i = 0; i < obj->keys.size(); ++i) {
        JsonObject *mat_spec = obj->values[i]->to_object();
        if (! mat_spec) {
            Log::error() << "Invalid material: " << obj->keys[i] << "." << endn;
            continue;
        }
                
        string mat_type = mat_spec->get_value("type")->to_string();
        if (mat_type == "lambertian") {
            JsonArray *a = mat_spec->get_value("kd")->to_array();
            if (!a) {
                Log::error() << "Invalid Lambertian kd: " << obj->keys[i] << "." << endn;
                continue;
            }

            Rgb kd = a->to_rgb();            
            Material *m = new BsdfMaterial(new Lambertian(kd));
            m->set_name(obj->keys[i]);
            materials->push_back(m);

        } else if (mat_type == "ward") {
            JsonArray *a_kd = mat_spec->get_value("kd")->to_array();
            if (!a_kd) {
                Log::error() << "Invalid Ward BRDF kd: " << obj->keys[i] << "." << endn;
                continue;
            }
            JsonArray *a_ks = mat_spec->get_value("ks")->to_array();
            if (!a_ks) {
                Log::error() << "Invalid Ward BRDF ks: " << obj->keys[i] << "." << endn;
                continue;
            }
            Float ax = mat_spec->get_value("ax")->to_numeric();
            Float ay = mat_spec->get_value("ay")->to_numeric();
            Rgb kd = a_kd->to_rgb();
            Rgb ks = a_ks->to_rgb();

            Material *m = new BsdfMaterial(new Ward(kd, ks, ax, ay));
            m->set_name(obj->keys[i]);
            materials->push_back(m);

        } else if (mat_type == "phong") {
            JsonArray *a_kd = mat_spec->get_value("kd")->to_array();
            if (!a_kd) {
                Log::error() << "Invalid Phong BRDF kd: " << obj->keys[i] << "." << endn;
                continue;
            }
            JsonArray *a_ks = mat_spec->get_value("ks")->to_array();
            if (!a_ks) {
                Log::error() << "Invalid Phong BRDF ks: " << obj->keys[i] << "." << endn;
                continue;
            }
            Float n = mat_spec->get_value("n")->to_numeric();
            Rgb kd = a_kd->to_rgb();
            Rgb ks = a_ks->to_rgb();

            Material *m = new BsdfMaterial(new ModifiedPhong(kd, ks, n));
            m->set_name(obj->keys[i]);
            materials->push_back(m);

        } else if (mat_type == "glass") {

            Material *m;
            Float n = mat_spec->get_value_strict("n")->to_numeric();

            if (mat_spec->has_key("kd")) {
                Rgb kd = mat_spec->get_value_strict("kd")->to_array()->to_rgb();            
                m = new BsdfMaterial(new Glass(n, kd));
            } else {
                m = new BsdfMaterial(new Glass(n));
            }
            m->set_name(obj->keys[i]);
            materials->push_back(m);

        } else if (mat_type == "mirror") {

            Material *m;
            if (mat_spec->has_key("kd")) {
                Rgb kd = mat_spec->get_value("kd")->to_array()->to_rgb();                
                m = new BsdfMaterial(new Mirror(kd));
            } else {
                m = new BsdfMaterial(new Mirror());
            }            
            m->set_name(obj->keys[i]);
            materials->push_back(m);

        } else {
            Log::error() << "Unknown material type: " << obj->keys[i] << "." << endn;
            continue;
        }
    }
}

void JsonImporter::load_light_specs(JsonObject *obj) {
    for (int i = 0; i < obj->keys.size(); ++i) {
        JsonObject *light_spec = obj->values[i]->to_object();
        if (! light_spec) {
            Log::error() << "Invalid light spec: " << obj->keys[i] << "." << endn;
            continue;
        }
                
        string type = light_spec->get_value("type")->to_string();
        if (type == "area") {
            JsonArray *a_rad = light_spec->get_value("radiance")->to_array();
            if (!a_rad) {
                Log::error() << "Invalid radiance in light spec: " << obj->keys[i] << "." << endn;
                continue;
            } else {
                Rgb radiance = a_rad->to_rgb();
                AreaLightSpec *spec = new AreaLightSpec;
                spec->set_name(obj->keys[i]);
                spec->type = Light::AREA_LIGHT;
                spec->radiance = radiance;
                light_specs->push_back(spec);
            }
        } else if (type == "directional") {
            Rgb radiance = light_spec->get_value_strict("radiance")->to_array()->to_vec3();
            Vec3 dir = light_spec->get_value_strict("direction")->to_array()->to_vec3();

            DirectionalLightSpec *spec = new DirectionalLightSpec;
            spec->set_name(obj->keys[i]);
            spec->type = Light::DIRECTIONAL_LIGHT;
            spec->radiance = radiance;
            spec->dir = dir;
            light_specs->push_back(spec);            

        } else if (type == "point") {
            JsonArray *a_rad = light_spec->get_value("radiance")->to_array();
            if (!a_rad) {
                Log::error() << "Invalid radiance in light spec: " << obj->keys[i] << "." << endn;
                continue;
            } 
            Rgb radiance = a_rad->to_rgb();

            JsonArray *a_pos = light_spec->get_value("position")->to_array();
            if (!a_pos) {
                Log::error() << "Invalid position in light spec: " << obj->keys[i] << "." << endn;
                continue;
            }
            Vec3 pos = a_pos->to_vec3();

            PointLightSpec *spec = new PointLightSpec;
            spec->set_name(obj->keys[i]);
            spec->type = Light::POINT_LIGHT;
            spec->pos = pos;
            spec->radiance = radiance;
            light_specs->push_back(spec);
        
        } else if (type == "spot") {
            JsonArray *a_rad = light_spec->get_value("radiance")->to_array();
            if (!a_rad) {
                Log::error() << "Invalid radiance in light spec: " << obj->keys[i] << "." << endn;
                continue;
            } 
            Rgb radiance = a_rad->to_rgb();

            JsonArray *a_pos = light_spec->get_value("position")->to_array();
            if (!a_pos) {
                Log::error() << "Invalid position in light spec: " << obj->keys[i] << "." << endn;
                continue;
            }
            Vec3 pos = a_pos->to_vec3();
            Vec3 dir = light_spec->get_value("direction")->to_array()->to_vec3();
            Float phi = light_spec->get_value("phi")->to_numeric();
            Float theta = light_spec->get_value("theta")->to_numeric();
            Float falloff = light_spec->get_value("falloff")->to_numeric();

            SpotLightSpec *spec = new SpotLightSpec;
            spec->set_name(obj->keys[i]);
            spec->type = Light::SPOT_LIGHT;
            spec->pos = pos;
            spec->dir = dir;
            spec->phi = phi;
            spec->theta = theta;
            spec->falloff = falloff;
            spec->radiance = radiance;
            light_specs->push_back(spec);

        } else {
            Log::error() << "Unknown light type: " << obj->keys[i] << "." << endn;
        }
    }
}

void JsonImporter::load_cameras(JsonObject *obj) {
    for (int i = 0; i < obj->keys.size(); ++i) {
        JsonObject *camera_spec = obj->values[i]->to_object();
        if (! camera_spec) {
            Log::error() << "Invalid camera: " << obj->keys[i] << "." << endn;
            continue;
        }

        string type = camera_spec->get_value("type")->to_string();
        if (type == "pinhole") {
            Vec3 pos = camera_spec->get_value("pos")->to_array()->to_vec3();
            Vec3 lookat = camera_spec->get_value("lookat")->to_array()->to_vec3();
            Vec3 up = camera_spec->get_value("up")->to_array()->to_vec3();
            Float near_plane = camera_spec->get_value("near")->to_numeric();
            Float far_plane = camera_spec->get_value("far")->to_numeric();
            Size2 img_size = camera_spec->get_value("image")->to_array()->to_size2();
            Size2 film_size = camera_spec->get_value("film")->to_array()->to_size2();

            JsonValue *val_focal = camera_spec->get_value("focal");
            Float focal;
            if (val_focal) {
                focal = val_focal->to_numeric();
            } else {
                JsonValue *val_fov = camera_spec->get_value("fov");
                if (!val_fov) {
                    Log::info() << "Focal/field of view parameter missing." << endn;
                } else {
                    Float fov = val_fov->to_numeric();
                    focal = Camera::focal_length_from_fov(fov, film_size.height);
                }
            }            

            Camera *camera = new Camera();
            Sensor *sensor = new Sensor(img_size, film_size);
            camera->set_name(obj->keys[i]);
            camera->set_perspective(pos, lookat, up, focal, sensor);            
            camera->set_near_plane(near_plane);
            camera->set_far_plane(far_plane);
            cameras->push_back(camera);
        } else {
            Log::error() << "Unknown camera type: " << obj->keys[i] << "." << endn;
        }
    }
}

void JsonImporter::load_scene(JsonObject *obj) {
    Log::info() << "Loading world ..." << endn;    
    JsonObject *world = obj->get_value_strict("world")->to_object();
    if (!world) {
        Log::error() << "Invalid world definition." << endn;
        return;
    }
    load_world(world);

    Log::info() << "Loading views ..." << endn;
    JsonObject *views = obj->get_value_strict("views")->to_object();
    if (!views) { 
        Log::error() << "Invalid view list." << endn;
        return;
    }
    load_camera_views(views);

    Log::info() << "Loading other options ..." << endn;
    JsonValue *scene_options_val = obj->get_value("options");
    if (scene_options_val != NULL) {
        scene_options = scene_options_val->to_object();
    } else {
        Log::info() << "No scene options found." << endn;
    }
}

void JsonImporter::load_world(JsonObject *obj) {
    Log::info() << "Loading world surfaces ..." << endn;
    JsonObject *surfaces = obj->get_value("surfaces")->to_object();
    if (!surfaces) {
        Log::error() << "Invalid world surfaces." << endn;
        return;
    }
    load_surfaces(surfaces);

    Log::info() << "Loading world lights ..." << endn;
    JsonObject *lights = obj->get_value("lights")->to_object();
    if (!lights) {
        Log::error() << "Invalid world lights." << endn;
        return;
    }
    load_lights(lights);
}

void JsonImporter::load_camera_views(JsonObject *views) {
    for (int i = 0; i < views->keys.size(); ++i) {
        JsonObject *view = views->values[i]->to_object();
        if (!view) {
            Log::error() << "Invalid view." << endn;
            continue;
        }
        
        JsonObject *integrator_spec = view->get_value("integrator")->to_object();
        if (!integrator_spec) { 
            Log::error() << "Invalid integrator." << endn;
            continue;
        }

        // to ease experiments, load integrator in scene file is deprecated. 
        Integrator *integrator = load_integrator(integrator_spec);
        if (!integrator) {
            continue;
        }

        string camera_name = view->get_value("camera")->to_string();

        Camera *camera = get_camera(camera_name);
        if (!camera) {
            Log::error() << "Camera " << camera_name << " not found." << endn;
            continue;
        }

        CameraView *camera_view = new CameraView();
        camera_view->set_integrator(integrator);
        camera_view->set_camera(camera);
        camera_view->set_name(view->keys[i]);
        camera_views->push_back(camera_view);
    }
}

Integrator *JsonImporter::load_integrator(JsonObject *spec) {
    int max_bounce = spec->get_value("max_bounce")->to_int();    
    if (max_bounce == 0) {
        Log::warn() << "Integrator max bounce is zero." << endn;
    }
    
    string type = spec->get_value("type")->to_string();
    
    if (type == "path" || type == "vpl") {
        int samples = spec->get_value("samples")->to_int();    
        if (samples == 0) {
            Log::warn() << "Monte Carlo integrator total sample is zero." << endn;
        }
        MonteCarloIntegrator *integrator = NULL;
        
        if (type == "path") 
            integrator = new PathTracing();
        else if (type == "vpl")
            integrator = new VirtualPointLight();        
        
        if (integrator == NULL)
            integrator = new PathTracing();
        
        integrator->set_max_bounce(max_bounce);
        integrator->set_num_samples(samples);

        return integrator;

    }
    /* 
    else if (type == "vpl_gpu") {
        
        int particles = spec->get_value("particles")->to_int();    
        if (particles == 0) {
            Log::warn() << "Virtual point light integrator total starting particle is zero." << endn;
        }

        bool direct = spec->get_value("direct")->to_bool();

        VirtualSphericalLightGpu *integrator = new VirtualSphericalLightGpu();
        integrator->set_max_bounce(max_bounce);    
        integrator->set_start_particles(particles);
        integrator->set_direct_lighting(direct);
        return integrator;
        
    } else if (type == "vsl_gpu") {
        /*
        int particles = spec->get_value("particles")->to_int();    
        if (particles == 0) {
            Log::warn() << "Virtual spherical light integrator total starting particle is zero." << endn;
        }

        bool direct = spec->get_value("direct")->to_bool();

        VirtualSphericalLightGpu *integrator = new VirtualSphericalLightGpu();
        integrator->set_max_bounce(max_bounce);    
        integrator->set_start_particles(particles);
        integrator->set_direct_lighting(direct);
        return integrator;
       
    }*/
    else {
        Log::error() << "Invalid integrator type: " << type << "." << endn;
        return NULL;
    }
}

void JsonImporter::load_surfaces(JsonObject *obj) {
    for (int i = 0; i < obj->keys.size(); ++i) {
        
        string name = obj->keys[i];

        if (obj->values[i]->is_string()) {

            // surfaces are defined with geometry and material in some scene files such as OBJ. 
            // just load them all
            string file = obj->values[i]->to_string();
            string obj_file = scene_path + "mesh/" + file;
            LogStream &out = Log::info();
            out << "Importing geometry and material from OBJ: " << setw(30) << left << file << " ... " << endn;
            Stats stats;
            stats.tic();
            ObjImporter obj_importer;
            obj_importer.set_home(scene_path);
            obj_importer.load(obj_file);
            Shapes local_shapes;            
            load_shapes_from_obj(name, obj_importer, local_shapes);
            Materials local_materials;
            obj_importer.get_materials(local_materials);
            Textures local_textures;
            obj_importer.get_textures(local_textures);
            Surfaces local_surfaces;
            obj_importer.get_surfaces(local_surfaces);
            Lights local_lights;                            // include all area lights if any
            obj_importer.get_lights(local_lights);

            // no names are assigned for these mesh and materials
            shapes->insert(    shapes->end(),    local_shapes.begin(),    local_shapes.end());
            materials->insert( materials->end(), local_materials.begin(), local_materials.end());
            textures->insert(  textures->end(),  local_textures.begin(),  local_textures.end());
            surfaces->insert(  surfaces->end(),  local_surfaces.begin(),  local_surfaces.end());
            lights->insert(    lights->end(),    local_lights.begin(),    local_lights.end());

            stats.toc();
            out << setw(10) << " ... " << stats.elapsed() << " ms" << endn;            
            
        } else {

            JsonObject *surf = obj->values[i]->to_object();            
            if (!surf) {
                Log::error() << "Invalid surface: " << name << "." << endn;
                continue;
            }

            if (surf->has_key("visible")) {
                bool visible = surf->get_value("visible")->to_bool();
                if (! visible) continue;
            }

            string mesh_name = surf->get_value("mesh")->to_string();
            if (mesh_name == "") {
                Log::error() << "Invalid geometry name in surface: " << name << "." << endn;
                continue;
            }
                    
            string brdf_name = surf->get_value("bsdf")->to_string();
            if (brdf_name == "") {
                Log::error() << "Invalid material name in surface: " << name << "." << endn;
                continue;
            }

            Shapes pieces;
            get_shape(mesh_name, pieces);
        
            Material *m = get_material(brdf_name);
            if (! m) {
                Log::error() << "Invalid material reference in surface: " << name << "." << endn;
                continue;
            }

            for (int j = 0; j < pieces.size(); ++j) {
                int k = pieces[j]->get_surface_index();     // if appended it should have a positive index
                if (k >= 0) {
                    
                    if ((*surfaces)[k].get_shape() == pieces[j]) 
                        (*surfaces)[k] = Surface(pieces[j], m);
                    else
                        Log::info() << "Surface/shape index error." << endn;

                } else {
                    Surface surface(pieces[j], m);
                    surfaces->append(surface);                    
                }
            }

        }
    }
}

void JsonImporter::load_lights(JsonObject *obj) {
    for (int i = 0; i < obj->keys.size(); ++i) {
        JsonValue *val = obj->values[i];
        if (val->is_string()) {
            string light_name = val->to_string();

            LightSpec *light_spec = get_light_spec(light_name);
            if (! light_spec) {
                Log::error() << "Invalid light reference: " << obj->keys[i] << " : " << light_name << "." << endn;
                continue;
            }

            if (light_spec->type == Light::AREA_LIGHT) {

                Log::error() << "Invalid area light reference. Mesh property required: " << obj->keys[i] << " : " << light_name << endn;
                continue;

            } else if (light_spec->type == Light::POINT_LIGHT) {

                PointLightSpec *spec = (PointLightSpec *)light_spec;
                lights->push_back(new PointLight(spec->pos, spec->radiance));

            } else if (light_spec->type == Light::SPOT_LIGHT) {

                SpotLightSpec *spec = (SpotLightSpec *)light_spec;
                SpotLight *spot_light = new SpotLight(spec->pos, 
                                                      unit_vector(spec->dir),
                                                      spec->phi * A_PI / 180.0f, 
                                                      spec->theta * A_PI / 180.0f, 
                                                      spec->falloff, 
                                                      spec->radiance);
                lights->push_back(spot_light);

            } else if (light_spec->type == Light::DIRECTIONAL_LIGHT) {

                DirectionalLightSpec *spec = (DirectionalLightSpec *)light_spec;
                DirectionalLight *dir_light = new DirectionalLight(
                                                        unit_vector(spec->dir),
                                                        spec->radiance);
                lights->push_back(dir_light);

            } else {

                Log::error() << "Unsupported light type: " << obj->keys[i] << " : " << light_name << "." << endn;
                continue;

            }

        } else if (val->is_object()) {                      // for area light only
            
            JsonObject *spec = obj->values[i]->to_object();
            if (! spec) {
                Log::error() << "Invalid light spec at: " << obj->keys[i] << endn;
                continue;
            }
            string name = obj->keys[i];
            string light_name = spec->get_value_strict("light")->to_string();
                
            LightSpec *light_spec = get_light_spec(light_name);
            if (! light_spec) {
                Log::error() << "Invalid light reference: " << name << " : " << light_name << endn;
                continue;
            }

            if (light_spec->type == Light::AREA_LIGHT) {

                string mesh_name = spec->get_value_strict("mesh")->to_string();            
                Shapes pieces;
                get_shape(mesh_name, pieces);

                for (int j = 0; j < pieces.size(); ++j) {                
                    Surface surface(pieces[j], light_spec->radiance);
                    surfaces->append(surface);
                    lights->push_back(surface.get_area_light());
                }

            } else {

                Log::error() << "Only area light can be bound to mesh: " << name << " : " << light_name << endn;
                continue;

            }
        }
    }
}

void JsonImporter::apply_scene_options(Scene *scene, JsonObject *options) {
    if (options == NULL) return;

    if (options->get_value("geometry")) {
        
        JsonObject *geometry_options = scene_options->get_value("geometry")->to_object();
        if (geometry_options->has_key("tmin")) {
            Float tmin = geometry_options->get_value("tmin")->to_numeric();
            scene->set_tmin(tmin);
        }

        if (geometry_options->has_key("two_sided_surfaces")) {            
            scene->set_two_sided_surfaces(geometry_options->get_value("two_sided_surfaces")->to_bool());
        }

    }


    JsonValue *opengl_val;
    if ((opengl_val = options->get_value("opengl")) != NULL) {
        JsonObject *opengl = opengl_val->to_object();

        JsonValue *val;
        if ((val = opengl->get_value("shadow_map")) != NULL) {
            JsonObject *sm_options = val->to_object();
            ShadowMapOptions shadow_map_options;
            shadow_map_options.polygon_offset_factor = sm_options->get_numeric("polygon_offset_factor");
            shadow_map_options.polygon_offset_unit   = sm_options->get_numeric("polygon_offset_unit");
            scene->set_shadow_map_options(shadow_map_options);
        }

        if ((val = opengl->get_value("sensor")) != NULL) {
            JsonObject *s_options = val->to_object();
            SensorOptions sensor_options;
            string s = s_options->get_string("multi_sample");
            if (s == "msaa")
                sensor_options.multi_sample = SensorOptions::MSAA;
            else
                sensor_options.multi_sample = SensorOptions::NONE;
            sensor_options.num_antialiased_samples = s_options->get_int("num_antialiased_samples");
            scene->set_sensor_options(sensor_options);
        }
    }

    JsonValue *experimental_val;
    if ((experimental_val = options->get_value("experimental")) != NULL) {
        JsonObject *experimental = experimental_val->to_object();

        ParamDict *dict = scene->get_param_dict();
        for (int i = 0; i < experimental->keys.size(); ++i) {
            JsonValue *val = experimental->values[i];
            if (val->type == JsonValue::JSON_VALUE_PRIMITIVE) {
                JsonPrimitive *prim = val->data.primitive;
                if (prim->is_float())
                    dict->add(experimental->keys[i], prim->to_float());
                if (prim->is_int())
                    dict->add(experimental->keys[i], prim->to_int());
                if (prim->is_string())
                    dict->add(experimental->keys[i], prim->to_string());                
            } else {
                Log::warn() << "Unsupported array/object for option: " << experimental->keys[i];
            }
        }
    }
}

void JsonImporter::load(const string path) {
    // entry point
    shapes = new Shapes();
    materials = new Materials();
    textures = new Textures();
    surfaces = new Surfaces();
    lights = new Lights();
    light_specs = new LightSpecs();
    cameras = new Cameras();
    camera_views = new CameraViews();

    string file_name;
    File::split_path(path, file_name, scene_path, scene_name);
	Log::info() << "Scene name: " << scene_name << endn;
	Log::info() << "Scene path: " << scene_path << endn;

    JsonParser *scene_parser = new JsonParser(this);
    scene_parser->load(path);
}

Scene *JsonImporter::get_scene() {
    Random *rd = new Random();
    
    Aggregate *agg = NULL;
    if (scene_options->get_value("geometry")) {
        
        JsonObject *geometry_options = scene_options->get_value("geometry")->to_object();
        if (geometry_options->has_key("traversal")) {

            string traversal_type = geometry_options->get_value("traversal")->to_string();
            if (traversal_type == "bvh") {

                agg = new AggregateBvh(*surfaces);
                
            } else if (traversal_type == "none") { 

                agg = new Aggregate(*surfaces);

            } else {    
                                
                Log::error() << "Invalid traversal data structure. Revert to BVH." << endn;
                agg = new AggregateBvh(*surfaces);
            }        
        }
        
    } 
    
    if (! agg) {    // default is BVH
        agg = new AggregateBvh(*surfaces);
    }

    CameraView *view = (*camera_views)[0];
    Camera *camera = view->get_camera();
    scene = new Scene(scene_name, agg, shapes, 
                      materials, textures, 
                      surfaces, lights, camera, rd);

    Integrator *integrator = view->get_integrator();
    scene->set_integrator(integrator);        

    apply_scene_options(scene, scene_options);

    return scene;
}

} // end namespace
