#ifndef _RAY_H_
#define _RAY_H_

#include "vec3.h"

#ifdef RZ_MORE_CHECK
#include "log.h"
#endif

using namespace std;

namespace Renzoku {

class Ray {
public:
    Ray() {}
    Ray(const Vec3 &_o, const Vec3 &_d) : o(_o), d(unit_vector(_d)) { 
        // avoid NaN as this will cause traversal performance dropped a few hundred times 
        if (d.is_nan())
            d = Vec3(0.0f, 0.0f, 0.0f);

        // avoid NaN here otherwise box hit checking is not accurate when one of elements in d is zero.                
        Vec3 dd(
            fabs(d.x()) < ZERO_EPSILON ? ZERO_EPSILON : d.x(),
            fabs(d.y()) < ZERO_EPSILON ? ZERO_EPSILON : d.y(),
            fabs(d.z()) < ZERO_EPSILON ? ZERO_EPSILON : d.z());
            
#ifdef RZ_MORE_CHECK
        if (d.is_nan() || (d.x() == 0.0f && d.y() == 0.0f && d.z() == 0.0f)) {
            Log::info() << "NaN/zero dir ray can cause ray traversal to be extremely slow." << endn;
        }
#endif

        inv_d = Vec3(1.0f) / dd; 
        org_inv_d = -o * inv_d; 
    }
    Ray(const Ray &r) { *this = r; }
    
    Vec3 org() const { return o; }
    Vec3 dir() const { return d; }
    Vec3 point(Float t) const { return o + t * d; }    

    friend ostream& operator<<(ostream &os, const Ray &r);

public:
    Vec3 o, d, inv_d, org_inv_d;
    Float n;        // the medium index the ray is travelling in
};

inline ostream& operator<<(ostream &os, const Ray &r) {
    os << "(" << r.org() << ") + t(" << r.dir() << ")";
    return os;
}

struct RayBundle {
    vector<Ray> rays;
};
typedef vector<RayBundle> RayBundles;

};
#endif

