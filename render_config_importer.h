#ifndef _RENDER_CONFIG_IMPORTER_H_
#define _RENDER_CONFIG_IMPORTER_H_

#include "scene_importer.h"
#include "json_parser.h"

#include "integrator.h"

namespace Renzoku {

class RenderConfigImporter : public JsonElementVisitor,
                             public FileImporter, 
                             public IntegratorImporter {
    
public:
    RenderConfigImporter();

    virtual void load(const string file);
    virtual Integrator *get_integrator();

public:
    virtual bool process_object(int stack_level, string last_key, JsonObject *obj);
    virtual bool process_array(string last_key, JsonArray *array);
    virtual void process_root(JsonObject *root);

private:
    void load_integrator(JsonObject *spec);
    Integrator *integrator;

private:
    enum Keyword {
        KW_INTEGRATOR,
        KW_TOTAL
    };

    static string keywords[KW_TOTAL];
};

}  // end namespace

#endif