#ifndef _LIGHT_TREE_H_
#define _LIGHT_TREE_H_

#include "common.h"
#include "boundingbox.h"
#include "cone.h"
#include "brdf_point_light.h"
#include "log.h"

namespace Renzoku {

struct LightNode {
    LightNode *childs;      // pointer to the two childs. Representative child is always at index 0.
    int num_childs;
    
    BrdfPointLight *light;  // the representative light, chosen randomly based on power distribution.
    Float scale;            // scale the power of the representative to the amount of all lights in the cluster

    BoundingBox bb;         // bounding box
    Cone bc;                // bounding cone

    int index;
    bool clustered;

    LightNode() : childs(NULL), num_childs(0), light(NULL), scale(1.0f), clustered(false), 
                  bc(Cone(Vec3(0.0f), 0.0f)) { // degenerated cone
    }

    ~LightNode() {
        
    }

    void find_bounding() {
        if (num_childs <= 0) {      // leaf node's bounding box
            bb = BoundingBox(light->org());
            bc = Cone(light->normal(), 0.0f);
            if (light->normal().is_nan()) {
                Log::info() << "Light normal NaN" << endn;
            }
            return;
        }

        bb = BoundingBox(childs[0].bb);
        if (num_childs > 1)
            bb.merge(childs[1].bb);

        bc = Cone(childs[0].bc);
        if (num_childs > 1)
            bc.merge_direction(childs[1].bc);

        if (bc.normal().is_nan()) {
            Log::info() << "Cone normal NaN" << endn;

            Cone bc2 = Cone(childs[0].bc);
            if (num_childs > 1)
                bc2.merge_direction(childs[1].bc);
        }
    }
};

}

#endif