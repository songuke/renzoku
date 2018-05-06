#ifndef _CONE_H_
#define _CONE_H_

#include "common.h"
#include "vec3.h"
#include "onb.h"

#include <cmath>
using namespace std;

namespace Renzoku {

class Cone {
public:	
    Cone(const Cone &c);
	Cone(const Vec3 &p, const Vec3 &n, Float half_angle);
	Cone(const Vec3 &n, Float half_angle);

	/**
	 * A cone spans from src to dest in anti-clockwise direction.
	 */
	Cone(const Vec3 &src, const Vec3 &dest);

	inline Vec3 org() const;
	inline Vec3 normal() const;
	inline Float half_angle() const;
	
public:
	/**
	 * Find the normal and half angle of the bounding cone that swipes from b to a in anti-clockwise order.
	 * 
	 * The position of the new cone is set to the first cone.
	 */
	static Cone merge_direction(const Cone &src, const Cone &dest);
        
    /**
     * The cone rotates from a source place to the current cone position.
     */
    void merge_direction(const Cone &src);
    void merge_direction_algebraic(const Vec3 &v);
    void merge_direction(const Vec3 &v);
    
    bool contains(const Vec3 &dir) const;

    inline void set_half_angle(Float half_angle);

private:
    /**
     * Angle to rotate src to dest in range [0, 360). 
     * 
     * Assume: src, dest are in XY plane.
     */
    static Float angle(const Vec3 &src, const Vec3 &dest);

    /**
     * Return the angle in range [0, 360) from x-axis to a vector in XY plane.
     */ 
    static Float angle_from_x_axis(const Vec3 &a);


private:
	Vec3 p;
	Vec3 n;
	Float half;
};

inline Cone::Cone(const Cone &c) : p(c.p), n(c.n), half(c.half) {
}

inline Cone::Cone(const Vec3 &p, const Vec3 &n, Float half_angle) : p(p), n(n), half(half_angle) {
}

inline Cone::Cone(const Vec3 &n, Float half_angle) : p(Vec3(0.f)), n(n), half(half_angle) {
}

inline Cone::Cone(const Vec3 &src, const Vec3 &dest) : p(Vec3(0.f)) {    
    Vec3 l = unit_vector(src);
    Vec3 r = unit_vector(dest);
    Onb uvn;
    uvn.init_from_uv(l, r); // anti-clockwise
    l = uvn.world_to_local(l);
    r = uvn.world_to_local(r);    

    half = angle(l, r) * 0.5f;

    n = uvn.local_to_world(unit_vector((l + r) / 2));    
    if (half > A_PI)
        n = -n;
}

inline Vec3 Cone::org() const {
	return p;
}

inline Vec3 Cone::normal() const {
	return n;
}

inline Float Cone::half_angle() const {
	return half;
}

inline void Cone::set_half_angle(Float half_angle) {
    half = half_angle;
}

inline bool Cone::contains(const Vec3 &v) const {
    return dot(v, n) >= cos(half) ? true : false;
}

} // end namespace

#endif