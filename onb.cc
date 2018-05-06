#include "onb.h"

namespace Renzoku {

void Onb::set_identity() {
    uu = Vec3(1.0f, 0.0f, 0.0f);
    vv = Vec3(0.0f, 1.0f, 0.0f);
    nn = Vec3(0.0f, 0.0f, 1.0f);
}

void Onb::init_from_u(const Vec3 &_u) {
    // vector vv and nn are arbitrary, so does not matter to cross with i or j.
    Vec3 i(1.0f, 0.0f, 0.0f);
    Vec3 j(0.0f, 1.0f, 0.0f);

    uu = unit_vector(_u);
    vv = cross(uu, i);
    if (vv.length() < ONB_EPSILON)
        vv = cross(uu, j);
    nn = cross(uu, vv);
    
    // make sure that vv is unit since uu and i or j are not orthogonal
    vv.make_unit_vector();
    nn.make_unit_vector();
}

void Onb::init_from_v(const Vec3 &_v) {
    Vec3 i(1.0f, 0.0f, 0.0f);
    Vec3 j(0.0f, 1.0f, 0.0f);

    vv = unit_vector(_v);
    uu = cross(vv, i);
    if (uu.length() < ONB_EPSILON)
        uu = cross(vv, j);
    nn = cross(uu, vv);

    // make sure that vv is unit since uu and i or j are not orthogonal
    uu.make_unit_vector();
    nn.make_unit_vector();
}

void Onb::init_from_n(const Vec3 &_n) {
    Vec3 i(1.0f, 0.0f, 0.0f);
    Vec3 j(0.0f, 1.0f, 0.0f);

    nn = unit_vector(_n);
    uu = cross(nn, i);
    if (uu.length() < ONB_EPSILON)
        uu = cross(nn, j);
    vv = cross(nn, uu);

    // make sure that vv is unit since uu and i or j are not orthogonal
    uu.make_unit_vector();
    vv.make_unit_vector();
}

void Onb::init_from_uv(const Vec3 &_u, const Vec3 &_v) {
    Vec3 c = unit_vector(cross(_u, _v));
    if (c.is_nan()) {
        init_from_u(_u);
        return;
    }

    uu = unit_vector(_u);
    nn = c;
    vv = cross(nn, uu);
}

void Onb::init_from_vu(const Vec3 &_v, const Vec3 &_u) {
    Vec3 c = unit_vector(cross(_u, _v));
    if (c.is_nan()) {
        init_from_v(_v);
        return;
    }

    vv = unit_vector(_v);
    nn = c;
    uu = cross(vv, nn);
}

void Onb::init_from_un(const Vec3 &_u, const Vec3 &_n) {
    Vec3 c = unit_vector(cross(_n, _u));
    if (c.is_nan()) {
        init_from_u(_u);
        return;
    }

    uu = unit_vector(_u);
    vv = c;
    nn = cross(uu, vv);    
}

void Onb::init_from_nu(const Vec3 &_n, const Vec3 &_u) {
    Vec3 c = unit_vector(cross(_n, _u));
    if (c.is_nan()) {
        init_from_n(_n);
        return;
    }

    nn = unit_vector(_n);
    vv = c;
    uu = cross(vv, nn);
}

void Onb::init_from_vn(const Vec3 &_v, const Vec3 &_n) {
    Vec3 c = unit_vector(cross(_v, _n));
    if (c.is_nan()) {
        init_from_v(_v);
        return;
    }

    vv = unit_vector(_v);
    uu = c;
    nn = cross(uu, vv);
}

void Onb::init_from_nv(const Vec3 &_n, const Vec3 &_v) {
    Vec3 c = unit_vector(cross(_v, _n));
    if (c.is_nan()) {
        init_from_n(_n);
        return;
    }

    nn = unit_vector(_n);
    uu = c;
    vv = cross(nn, uu);
}

void Onb::flip_n() {
    nn = -nn;
    vv = -vv;    
}

Vec3 Onb::world_to_local(const Vec3 &w) const {
    return Vec3(dot(uu, w), dot(vv, w), dot(nn, w));
}

Vec3 Onb::local_to_world(const Vec3 &v) const {
    return v[0] * uu + v[1] * vv + v[2] * nn;
}

bool operator==(const Onb &o1, const Onb &o2) {
    return (o1.u() == o2.u() && o1.v() == o2.v() && o1.n() == o2.n()); }

istream& operator>>(istream &is, Onb &o) {
    Vec3 u, v, n;
    is >> u >> v >> n;
    o.init_from_uv(u, v);
    return is;
}

ostream& operator<<(ostream &os, const Onb &o) {
    os << o.u() << '\n' << o.v() << '\n' << o.n() << '\n';
    return os;
}

Vec3 Onb::spherical_to_world(Float theta, Float phi) const {
    Float sin_theta = sin(theta);
    Float cos_theta = cos(theta);
    Float sin_phi = sin(phi);
    Float cos_phi = cos(phi);
    Vec3 v(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta);
    return local_to_world(v);
}

void Onb::world_to_spherical(const Vec3 &v, Float &theta, Float &phi) const {
    Vec3 u = world_to_local(v);
    Float x = u.x();
    Float y = u.y();
    Float z = u.z();

    theta = acos(z);
    if (x == 0) {
        if (y >= 0) {
            phi = HALF_PI;
        } else {
            phi = A_PI + HALF_PI;
        }
    } else {
        phi = atan2(y, x);
        if (phi < 0) phi += 2 * A_PI; // shift to [0, 2 * M_PI]
    }
}

} // end namespace Renzoku
