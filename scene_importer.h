#ifndef _SCENE_IMPORTER_H_
#define _SCENE_IMPORTER_H_

#include "common.h"
#include "surface.h"

namespace Renzoku {

class FileImporter {
public:
    virtual void load(const string file) = 0;
};

class MaterialImporter {    
public:
    virtual void get_materials(Materials &materials) = 0;
};

class TextureImporter {
public:
    virtual void get_textures(Textures &textures) = 0;
};

class ShapeImporter {
public:
    virtual void get_shapes(Shapes &shapes) = 0;
};

class SurfaceImporter {
public:
    virtual void get_surfaces(Surfaces &surfaces) = 0;
};

class SceneImporter {
public:
    virtual Scene *get_scene() = 0;
};

class IntegratorImporter {
public:
    virtual Integrator *get_integrator() = 0;
};

} // end namespace

#endif