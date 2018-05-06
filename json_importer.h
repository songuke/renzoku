#ifndef _JSON_IMPORTER_H_
#define _JSON_IMPORTER_H_

#include "common.h"
#include "math3.h"
#include "named_object.h"

#include "json_parser.h"

#include "scene_importer.h"
#include "light.h"

namespace Renzoku {

struct LightSpec : NamedObject {
    Light::Type type;
    Rgb radiance;
};

struct AreaLightSpec : LightSpec {
};

struct PointLightSpec : LightSpec {
    Vec3 pos;
};

struct SpotLightSpec : LightSpec {
    Vec3 pos;
    Vec3 dir;
    Float phi, theta, falloff;
};

struct DirectionalLightSpec : LightSpec {
    Vec3 dir;
};

typedef vector<LightSpec *> LightSpecs;

class JsonImporter : public JsonElementVisitor,
                     public FileImporter, 
                     public SceneImporter {

public:
    JsonImporter();

    virtual void load(const string file);
    virtual Scene *get_scene();

public:    
    bool process_object(int stack_level, string last_key, JsonObject *obj);
    bool process_array(string last_key, JsonArray *array);
    void process_root(JsonObject *root);

private:
    void load_geometry(JsonObject *obj);
    void load_materials(JsonObject *obj);
    void load_light_specs(JsonObject *obj);
    void load_cameras(JsonObject *obj);
    void load_scene(JsonObject *obj);
    void load_surfaces(JsonObject *obj);
    void load_lights(JsonObject *obj);
    void load_world(JsonObject *obj);
    void load_camera_views(JsonObject *obj);
    Integrator *load_integrator(JsonObject *obj);
    void apply_scene_options(Scene *scene, JsonObject *options);

private:
    void get_shape(string name, Shapes &matched);
    Material *get_material(string name);
    LightSpec *get_light_spec(string name);
    Camera *get_camera(string name);
    
private:
    Shapes *shapes;
    Materials *materials;
    Textures *textures;
    LightSpecs *light_specs;
    Lights *lights;
    Cameras *cameras;
    CameraViews *camera_views;
    Surfaces *surfaces;
    JsonObject *scene_options;
    Scene *scene;
    string scene_name;
    string scene_path;

private:
    enum Keyword {
        KW_GEOMETRY,
        KW_MATERIALS,
        KW_LIGHT_SPECS,
        KW_CAMERAS,
        KW_SCENE,
        KW_WORLD,
        KW_SURFACES,
        KW_LIGHTS,
        KW_INTEGRATOR,
        KW_CAMERA,
        KW_TYPE,
        KW_TOTAL
    };
    static string keywords[KW_TOTAL];
};

} // end namespace

#endif
