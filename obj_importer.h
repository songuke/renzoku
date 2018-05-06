#ifndef _OBJ_IMPORTER_H_
#define _OBJ_IMPORTER_H_

#include "scene_importer.h"

#include "obj_parser.h"
#include "common.h"

namespace Renzoku {

/** 
 * Wavefront OBJ and MTL importer.
 * 
 * This is a wrapper class to Yingchao's customized OBJ loader, which
 * extends OBJ to store a complete scene (with light and camera added).
 */
class ObjImporter : public FileImporter, 
                    public ShapeImporter, 
                    public MaterialImporter, 
                    public TextureImporter,
                    public SurfaceImporter,
                    public SceneImporter {

public:
    ObjImporter();
    ~ObjImporter();

    virtual void load(const string file);

    /**
     * Standard OBJ file loader should only retrieve shapes, materials, and surfaces.
     */
    virtual void get_shapes(Shapes &shapes);
    virtual void get_materials(Materials &materials);
    virtual void get_textures(Textures &textures);
    virtual void get_surfaces(Surfaces &surfaces);
    
    /**
     * Non-standard OBJ file can store a scene.
     */
    virtual Scene *get_scene();

    /**
     * Our OBJ might carry some area lights.
     */
    virtual void get_lights(Lights &lights);

    /**
     * Object names can be referred in the JSON scene file. 
     */
    virtual void get_object_names(vector<string> &object_names);

    void set_home(const string home);

private:
    string scene_path;
    string name;
    Shapes *shapes;
    Materials *materials;
    Textures *textures;
    Surfaces *surfaces;
    Camera *camera;
    Lights *lights;
    Scene *scene;
    vector<string> object_names;

private:
	obj_scene_data obj_data;
};

} // end namespace

#endif