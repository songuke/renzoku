#ifndef _LIGHT_PARTICLE_H_
#define	_LIGHT_PARTICLE_H_

#include "photon.h"
#include "brdf_point_light.h"
#include "path.h"

namespace Renzoku {

struct DensityPoint {
    Vec3 p, n;
    Float count;

    DensityPoint(Vec3 p, Vec3 &n) : p(p), n(n), count(0.f) {
    }
};

typedef vector<DensityPoint *> DensityPoints;

struct LightNode;

/*
 * An adapter class that encapsulates a photon or a virtual point light
 * so that kd-tree can treat both of them as light particles. 
 */    
class LightParticle {
public:
    LightParticle();
    LightParticle(Photon *photon);
    LightParticle(BrdfPointLight *light);
    LightParticle(DensityPoint *point);
    LightParticle(PathNode *node);
    LightParticle(LightNode *node);

public:
    inline Vec3 org() const;
    
    inline Photon *get_photon() const;
    inline BrdfPointLight *get_brdf_point_light() const;
    inline DensityPoint *get_density_point() const;
    inline PathNode *get_path_node() const;
    inline LightNode *get_light_node() const;

protected:
    union {
        Photon *photon;
        BrdfPointLight *light;
        DensityPoint *point;
        PathNode *node;                    // path node in path tracing
        LightNode *light_node;             // light tree node in Lightcuts
    };

    Vec3 p;    
};

typedef vector<LightParticle>         LightParticles;
typedef vector<LightParticle *>       LightParticlePtrs;

/**
 * This max heap is used to store the results of nearest neighbor query from a point p.
 * 
 * The light particles are sorted from furthest 
 * to nearest distance between p and the particle location.
 *
 * Since p is different for each query, we need to initialize 
 * the priority queue with a functor that takes p as a constructor parameter.
 */
struct LightParticleLess {
    LightParticleLess(const Vec3 &p) : p(p) {}

    bool operator()(LightParticle *a, LightParticle *b) {
        return (a->org() - p).squared_length() < (b->org() - p).squared_length();
    }

    Vec3 p;     // the query point
};

struct LightParticleHeap 
        : priority_queue<LightParticle*, LightParticlePtrs, LightParticleLess> {

    LightParticleHeap(const Vec3 &p) 
        : priority_queue<LightParticle*, LightParticlePtrs, LightParticleLess>(LightParticleLess(p)) {
    }
};


/**
 * Functor to compare two "pointers" of LightParticle. 
 * Priority given to "larger" pointer (the top of the heap is greatest).
 */
/*
static bool operator<(const LightParticle &p1, const LightParticle &p2) {
    // We want to put the photon with furthest distance to p to the top of the heap 
    // (so we can quickly discard it if needed).
    if (p1.dist == -1 || p2.dist == -1) {
        cout << "Error photon comparison." << endl;
    }
    return p1.dist < p2.dist;
}
template <class T>
struct ptr_less : public binary_function<T, T, bool> {
    bool operator()(const T &a, const T &b) {
        return (*a) < (*b); // this will call operator>() function of type T
    }
};
typedef priority_queue<LightParticle*, LightParticles, ptr_less<LightParticle*> > LightParticleHeap;
*/


inline Vec3 LightParticle::org() const {
    return p;
}

inline Photon *LightParticle::get_photon() const {
    return photon;
}

inline BrdfPointLight *LightParticle::get_brdf_point_light() const {
    return light;
}

inline DensityPoint *LightParticle::get_density_point() const {
    return point;
}

inline PathNode *LightParticle::get_path_node() const {
    return node;
}

inline LightNode *LightParticle::get_light_node() const {
    return light_node;
}

} // end namespace

#endif

