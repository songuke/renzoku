#ifndef _HIT_RECORD_H_
#define _HIT_RECORD_H_

#include "common.h"
#include "vec3.h"
#include "onb.h"

namespace Renzoku {

/**
 * Minimal data for visibility test and nearest intersection tracking.
 */
class GeometryHit : public BaseObject {
public:    
    Float t; 
    bool hit;
    Shape *shape;
    Float beta, gamma;

    GeometryHit() : shape(NULL), hit(false), t(RAY_TMAX), beta(0.0f), gamma(0.0f) {             
    }
};

class HitRecord : public BaseObject {
public:
    bool hit;
    Float t, beta, gamma;

    Vec3 tangent, normal;
    Vec3 shading_normal;
    Onb uvn;
    
    // TODO: to use union here; we need a variable to indicate type
    Material *material;   
    AreaLight *light;
    EnvLight *env_light;    

    Shape *shape;
    Surface *surface;

    Vec2 uv;    // texture coordinates
        
    HitRecord() : material(NULL), light(NULL), env_light(NULL), shape(NULL), surface(NULL),
                  hit(false), t(RAY_TMAX), beta(0.0f), gamma(0.0f) {
             
    }

    /**
     * Copy minimal hit info for nearest hit point tracking.
     * Other info can be filled after the nearest hit point is found.
     */
    void copy_geometry_hit(const GeometryHit &r) {
        t = r.t;
        beta = r.beta;
        gamma = r.gamma;
        hit = r.hit;
        shape = r.shape;
    }

};

typedef vector<HitRecord> HitRecords;

} // end namespace Renzoku

#endif
