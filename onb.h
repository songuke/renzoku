#ifndef _ONB_H_
#define _ONB_H_

#include "base_object.h"
#include "vec3.h"

namespace Renzoku {

/**
 * Orthonormal basis. Right-handed. Use n as up vector.
 *
 * In XYZ coordinate system, Z is up vector. This is different from the camera where Y is often the up vector.
 */
class Onb : public BaseObject {
public:
    Onb() {}
    Onb(const Vec3 &u, const Vec3 &v, const Vec3 &n) : uu(u), vv(v), nn(n) {}
    
    void init_from_u(const Vec3 &_u);
    void init_from_v(const Vec3 &_v);
    void init_from_n(const Vec3 &_n);

    Vec3 u() const { return uu; }
    Vec3 v() const { return vv; }
    Vec3 n() const { return nn; }
    void set(const Vec3 &_u, const Vec3 &_v, const Vec3 &_n) { uu = _u; vv = _v; nn = _n; }
    void set_identity();

    void init_from_uv(const Vec3 &_u, const Vec3 &_v);
    void init_from_vu(const Vec3 &_v, const Vec3 &_u);
    
    void init_from_un(const Vec3 &_u, const Vec3 &_n);
    void init_from_nu(const Vec3 &_n, const Vec3 &_u);
        
    void init_from_vn(const Vec3 &_v, const Vec3 &_n);
    void init_from_nv(const Vec3 &_n, const Vec3 &_v);

    Vec3 world_to_local(const Vec3 &w) const;
    Vec3 local_to_world(const Vec3 &v) const;

    void flip_n();

    friend bool operator==(const Onb &o1, const Onb &o2);
    
    friend istream& operator>>(istream &is, Onb &t);
    friend ostream& operator<<(ostream &os, const Onb &t);
    
    /**
     * Move v to local frame and transform to spherical coordinates.
     *
     * N as up vector.
     */
    void world_to_spherical(const Vec3 &v, Float &theta, Float &phi) const;
    Vec3 spherical_to_world(Float theta, Float phi) const;

protected:
    Vec3 uu, vv, nn;
};

}  // end namespace

#endif

