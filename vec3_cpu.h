#ifndef _VEC3_X87_H_
#define _VEC3_X87_H_

#include "common.h"
#include "random.h"

#include <math.h>
#include <iostream>
using namespace std;

namespace Renzoku {
class Vec3 {
public:
    Vec3();
    Vec3(Float e0, Float e1, Float e2);
    Vec3(Float v);
    Vec3(const Vec3 &v) { e[0] = v.e[0]; e[1] = v.e[1]; e[2] = v.e[2]; }
    
    Float x() const { return e[0]; }
    Float y() const { return e[1]; }
    Float z() const { return e[2]; }
    
    const Vec3& operator+() const;
    Vec3 operator-() const; 
    
    Float operator[](int i) const { return e[i]; }
    Float& operator[](int i) { return e[i]; }

    Float length() const;
    Float squared_length() const;
    
    void make_unit_vector();

    bool is_nan() const;
    friend bool is_nan(const Vec3 &v);

    void set_x(Float _x) { e[0] = _x; }
    void set_y(Float _y) { e[1] = _y; }
    void set_z(Float _z) { e[2] = _z; }

    Float min_component() const;
    Float max_component() const;
    Float min_abs_component() const;
    Float max_abs_component() const;

    int arg_min_component() const;
    int arg_max_component() const;
    int arg_min_abs_component() const;
    int arg_max_abs_component() const;
    
    friend istream& operator>>(istream &is, Vec3 &t);
    friend ostream& operator<<(ostream &os, const Vec3 &t);
  
    friend Vec3 operator+(const Vec3 &v1, const Vec3 &v2);
    //const Vec3 operator+(const Vec3 &v) const;            // any version (friend or member function) works the same.
    friend Vec3 operator-(const Vec3 &v1, const Vec3 &v2);
    friend Vec3 operator/(const Vec3 &v, Float t);
    friend Vec3 operator/(const Vec3 &u, const Vec3 &v);
    friend Vec3 operator*(const Vec3 &v, Float t);
    friend Vec3 operator*(Float scalar, const Vec3 &v);
    friend Vec3 operator*(const Vec3 &u, const Vec3 &v);
    
    Vec3& mul_add(const Vec3 &a, const Vec3 &b);
    friend Vec3 mul_add(const Vec3 &v, const Vec3& a, const Vec3 &b);

    Vec3& operator=(const Vec3 &v);
    Vec3& operator+=(const Vec3 &v);
    Vec3& operator-=(const Vec3 &v);
    Vec3& operator*=(const Float t);
    Vec3& operator/=(const Float t);

    friend bool operator<(const Vec3 &v, Float a);
    friend bool operator<(const Vec3 &v1, const Vec3 &v2);
    friend bool operator<=(const Vec3 &v, Float a);
    friend bool operator<=(const Vec3 &v1, const Vec3 &v2);
    friend bool operator==(const Vec3 &v, Float a);
    friend bool operator==(const Vec3 &v1, const Vec3 &v2);
    friend bool operator!=(const Vec3 &v, Float a);
    friend bool operator!=(const Vec3 &v1, const Vec3 &v2);
    friend bool operator>(const Vec3 &v, Float a);
    friend bool operator>(const Vec3 &v1, const Vec3 &v2);
    friend bool operator>=(const Vec3 &v, Float a);
    friend bool operator>=(const Vec3 &v1, const Vec3 &v2);
    
    bool operator<(const Vec3 &v);
    bool operator<(Float a);
    bool operator<=(const Vec3 &v);
    bool operator<=(Float a);
    bool operator==(const Vec3 &v);
    bool operator==(Float a);
    bool operator!=(const Vec3 &v);
    bool operator!=(Float a);
    bool operator>(const Vec3 &v);
    bool operator>(Float a);
    bool operator>=(const Vec3 &v);
    bool operator>=(Float a);
    
    friend Vec3 unit_vector(const Vec3 &v);
    friend Vec3 min(const Vec3 &v1, const Vec3 &v2);
    friend Vec3 max(const Vec3 &v1, const Vec3 &v2);
    friend Vec3 cross(const Vec3 &v1, const Vec3 &v2);
    friend Float dot(const Vec3 &v1, const Vec3 &v2);
    friend Float triple(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3);
    
    void clamp(const Vec3 &a, const Vec3 &b);
    friend Vec3 clamp(const Vec3 &u, const Vec3 &a, const Vec3 &b);

    /**
     * NOTE: this returns a vector in coordinate system where z-axis is up-vector.
     */
    static Vec3 from_spherical(Float theta, Float phi);
    void to_spherical(Float &theta, Float &phi);

    void random(Random &rd) { e[0] = rd(); e[1] = rd(); e[2] = rd(); }

    void rotate_about_x(Float angle);
    void rotate_about_y(Float angle);
    void rotate_about_z(Float angle);
    friend Vec3 rotate_about_x(const Vec3 &v, Float angle);
    friend Vec3 rotate_about_y(const Vec3 &v, Float angle);
    friend Vec3 rotate_about_z(const Vec3 &v, Float angle);

public:
    Float e[3];
    //Float e[4];       // NOTE: storing four components will make each frame render slower about 30ms due to copy overhead.
};

inline Vec3::Vec3() {
    e[0] = 0;
    e[1] = 0;
    e[2] = 0;
}

inline Vec3::Vec3(Float v) {
    //e[0] = e[1] = e[2] = v;
    e[0] = v;
    e[1] = v;
    e[2] = v;
}

inline Vec3::Vec3(Float e0, Float e1, Float e2) { 
    e[0] = e0; e[1] = e1; e[2] = e2; 
}

inline const Vec3& Vec3::operator+() const { return *this; }

inline Vec3 Vec3::operator-() const { return Vec3(-e[0], -e[1], -e[2]); }

inline Float Vec3::length() const { 
    return sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]); }

inline Float Vec3::squared_length() const {
    return e[0]*e[0] + e[1]*e[1] + e[2]*e[2]; }

inline void Vec3::make_unit_vector() {
    Float inv_len = 1.0f / length();
    e[0] *= inv_len;
    e[1] *= inv_len;
    e[2] *= inv_len;    
}

inline bool Vec3::is_nan() const {
    return ((e[0] != e[0]) || (e[1] != e[1]) || (e[2] != e[2]));
}

inline bool is_nan(const Vec3 &v) {
    return ((v.x() != v.x()) || (v.y() != v.y()) || (v.z() != v.z()));
}

inline Float Vec3::min_component() const {
    Float tmp = e[0];
    if (e[1] < tmp) tmp = e[1];
    if (e[2] < tmp) tmp = e[2];
    return tmp; 
}

inline Float Vec3::max_component() const {
    Float tmp = e[0];
    if (e[1] > tmp) tmp = e[1];
    if (e[2] > tmp) tmp = e[2];
    return tmp;
}

inline Float Vec3::min_abs_component() const {
    Float tmp = fabs(e[0]);
    if (fabs(e[1]) < tmp) tmp = fabs(e[1]);
    if (fabs(e[2]) < tmp) tmp = fabs(e[2]);
    return tmp;
}

inline Float Vec3::max_abs_component() const {
    Float tmp = fabs(e[0]);
    if (fabs(e[1]) > tmp) tmp = fabs(e[1]);
    if (fabs(e[2]) > tmp) tmp = fabs(e[2]);
    return tmp;
}

inline int Vec3::arg_min_component() const {
    int index = 0;
    Float tmp = e[0];
    if (tmp < e[1]) { tmp = e[1]; index = 1; }
    if (tmp < e[2]) { index = 2; }
    return index;
}

inline int Vec3::arg_max_component() const {
    int index = 0;
    Float tmp = e[0];
    if (tmp > e[1]) { tmp = e[1]; index = 1; }
    if (tmp > e[2]) { index = 2; }
    return index;
}

inline int Vec3::arg_min_abs_component() const {
    int index = 0;
    Float tmp = fabs(e[0]);
    if (tmp < fabs(e[1])) { tmp = fabs(e[1]); index = 1; }
    if (tmp < fabs(e[2])) { index = 2; }
    return index;
}

inline int Vec3::arg_max_abs_component() const {
    int index = 0;
    Float tmp = fabs(e[0]);
    if (tmp > fabs(e[1])) { tmp = fabs(e[1]); index = 1; }
    if (tmp > fabs(e[2])) { index = 2; }
    return index;
}

// friend implementation
inline Vec3 operator+(const Vec3 &v1, const Vec3 &v2) {
    return Vec3(v1.e[0] + v2.e[0], v1.e[1] + v2.e[1], v1.e[2] + v2.e[2]); 
}

// class function implementation
/*
inline const Vec3 Vec3::operator+(const Vec3 &v) const {
    return Vec3(e[0] + v.e[0], e[1] + v.e[1], e[2] + v.e[2]);
}*/

inline Vec3 operator-(const Vec3 &v1, const Vec3 &v2) {
    return Vec3(v1.e[0] - v2.e[0], v1.e[1] - v2.e[1], v1.e[2] - v2.e[2]);
}

inline Vec3 operator/(const Vec3 &v, Float t) {
    Float inv_t = 1.0f / t;
    return Vec3(v.e[0] * inv_t, v.e[1] * inv_t, v.e[2] * inv_t);
    //return Vec3(v.e[0]/t, v.e[1]/t, v.e[2]/t);        // Visual Studio 2012 does not automatically lift and compute 1/t first, so this line is ~0.5x slower.
}
    
inline Vec3 operator/(const Vec3 &u, const Vec3 &v) {
    return Vec3(u.e[0]/v.e[0], u.e[1]/v.e[1], u.e[2]/v.e[2]); }

inline Vec3 operator*(const Vec3 &v, Float t) {
    return Vec3(v.e[0]*t, v.e[1]*t, v.e[2]*t); }

inline Vec3 operator*(Float t, const Vec3 &v) {
    return Vec3(v.e[0]*t, v.e[1]*t, v.e[2]*t); }

inline Vec3 operator*(const Vec3 &u, const Vec3 &v) {
        return Vec3(u.e[0]*v.e[0], u.e[1]*v.e[1], u.e[2]*v.e[2]); }

inline Vec3& Vec3::operator=(const Vec3 &v) {
    e[0] = v.e[0]; e[1] = v.e[1]; e[2] = v.e[2]; return *this; }

inline Vec3& Vec3::operator+=(const Vec3 &v) {
    *this = *this + v; return *this; }

inline Vec3& Vec3::operator-=(const Vec3 &v) {
    *this = *this - v; return *this; }

inline Vec3& Vec3::operator*=(Float t) {
    *this = *this * t; return *this; }

inline Vec3& Vec3::operator/=(Float t) {
    *this = *this / t; return *this; }

inline Vec3& Vec3::mul_add(const Vec3& a, const Vec3& b) {
    e[0] = e[0] * a.e[0] + b.e[0];
    e[1] = e[1] * a.e[1] + b.e[1];
    e[2] = e[2] * a.e[2] + b.e[2];
    return *this;
}

inline Vec3 mul_add(const Vec3 &v, const Vec3& a, const Vec3 &b) {
    /*
    __m128 mv, ma, mb, mo;
    mv = _mm_load_ps(v.e);
    ma = _mm_load_ps(a.e);
    mb = _mm_load_ps(b.e);
    mo = _mm_fmadd_ps(mv, ma, mb);
    return Vec3(mo.m128_f32[1], mo.m128_f32[2], mo.m128_f32[3]);
    */
    
    return Vec3(v.e[0] * a.e[0] + b.e[0],
                v.e[1] * a.e[1] + b.e[1],
                v.e[2] * a.e[2] + b.e[2]);
}

inline Float dot(const Vec3 &v1, const Vec3 &v2) {
    return v1.e[0]*v2.e[0] + v1.e[1]*v2.e[1] + v1.e[2]*v2.e[2]; }

inline Vec3 cross(const Vec3 &v1, const Vec3 &v2) {
    return Vec3(
        v1.y() * v2.z() - v1.z() * v2.y(),
        v1.z() * v2.x() - v1.x() * v2.z(),
        v1.x() * v2.y() - v1.y() * v2.x()); 
}

inline bool operator<(const Vec3 &v, Float a) {
    return ((v.x() < a) && (v.y() < a) && (v.z() < a));
}

inline bool operator<(const Vec3 &v1, const Vec3 &v2) {
    return ((v1.x() < v2.x()) && (v1.y() < v2.y()) && (v1.z() < v2.z()));
}

inline bool operator<=(const Vec3 &v, Float a) {
    return ((v.x() <= a) && (v.y() <= a) && (v.z() <= a));
}

inline bool operator<=(const Vec3 &v1, const Vec3 &v2) {
    return ((v1.x() <= v2.x()) && (v1.y() <= v2.y()) && (v1.z() <= v2.z()));
}

inline bool operator==(const Vec3 &v, Float a) {
    return ((v.x() == a) && (v.y() == a) && (v.z() == a));
}

inline bool operator==(const Vec3 &v1, const Vec3 &v2) {
    return ((v1.x() == v2.x()) && (v1.y() == v2.y()) && (v1.z() == v2.z()));
}

inline bool operator!=(const Vec3 &v, Float a) {
    return ((v.x() != a) || (v.y() != a) || (v.z() != a));
}

inline bool operator!=(const Vec3 &v1, const Vec3 &v2) {
    return ((v1.x() != v2.x()) || (v1.y() != v2.y()) || (v1.z() != v2.z()));
}

inline bool operator>(const Vec3 &v, Float a) {
    return ((v.x() > a) && (v.y() > a) && (v.z() > a));
}
    
inline bool operator>(const Vec3 &v1, const Vec3 &v2) {
    return ((v1.x() > v2.x()) && (v1.y() > v2.y()) && (v1.z() > v2.z()));
}

inline bool operator>=(const Vec3 &v, Float a) {
    return ((v.x() >= a) && (v.y() >= a) && (v.z() >= a));
}

inline bool operator>=(const Vec3 &v1, const Vec3 &v2) {
    return ((v1.x() >= v2.x()) && (v1.y() >= v2.y()) && (v1.z() >= v2.z()));
}
    
inline bool Vec3::operator<(const Vec3 &v) {
    return ((e[0] < v.x()) && (e[1] < v.y()) && (e[2] < v.z()));
}
    
inline bool Vec3::operator<(Float a) {
    return ((e[0] < a) && (e[1] < a) && (e[2] < a));
}
    
inline bool Vec3::operator<=(const Vec3 &v) {
    return ((e[0] <= v.x()) && (e[1] <= v.y()) && (e[2] <= v.z()));
}
    
inline bool Vec3::operator<=(Float a) {
    return ((e[0] <= a) && (e[1] <= a) && (e[2] <= a));
}

inline bool Vec3::operator==(const Vec3 &v) {
    return ((e[0] == v.x()) && (e[1] == v.y()) && (e[2] == v.z()));
}
    
inline bool Vec3::operator==(Float a) {
    return ((e[0] == a) && (e[1] == a) && (e[2] == a));
}

inline bool Vec3::operator!=(const Vec3 &v) {
    return ((e[0] != v.x()) || (e[1] != v.y()) || (e[2] != v.z())); 
}

inline bool Vec3::operator!=(Float a) {
    return ((e[0] != a) || (e[1] != a) || (e[2] != a));
}
    
inline bool Vec3::operator>(const Vec3 &v) {
    return ((e[0] > v.x()) && (e[1] > v.y()) && (e[2] > v.z()));
}
    
inline bool Vec3::operator>(Float a) {
    return ((e[0] > a) && (e[1] > a) && (e[2] > a));
}
    
inline bool Vec3::operator>=(const Vec3 &v) {
    return ((e[0] >= v.x()) && (e[1] >= v.y()) && (e[2] >= v.z()));
}

inline bool Vec3::operator>=(Float a) {
    return ((e[0] >= a) && (e[1] >= a) && (e[2] >= a));
}

inline Vec3 unit_vector(const Vec3 &v) {    
    //return v / v.length();                  // this takes 3 divisions. 
    //Float inv_len = 1.0 / v.length();       // this is slower as the computation is in double data type.

    Float inv_len = 1.0f / v.length();
    //return v * inv_len;
    return Vec3(v.e[0] * inv_len, v.e[1] * inv_len, v.e[2] * inv_len);
}

inline Vec3 min(const Vec3 &v1, const Vec3 &v2) {
    Vec3 tmp(v1);
    if (v2.x() < v1.x()) tmp.set_x(v2.x());
    if (v2.y() < v1.y()) tmp.set_y(v2.y());
    if (v2.z() < v1.z()) tmp.set_z(v2.z());
    return tmp;
}

inline Vec3 max(const Vec3 &v1, const Vec3 &v2) {
    Vec3 tmp(v1);
    if (v2.x() > v1.x()) tmp.set_x(v2.x());
    if (v2.y() > v1.y()) tmp.set_y(v2.y());
    if (v2.z() > v1.z()) tmp.set_z(v2.z());
    return tmp;
}

inline Float triple(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3) {
    return dot(cross(v1, v2), v3); }

inline Vec3 from_spherical(Float theta, Float phi) {
    Float sin_theta = sin(theta);
    Float cos_theta = cos(theta);
    Float sin_phi = sin(phi);
    Float cos_phi = cos(phi);
    //return Vec3(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta);

    return Vec3(sin_theta * cos_phi, cos_theta, sin_theta * sin_phi); // y as up-vector
}

inline void Vec3::to_spherical(Float &theta, Float &phi) {    
    /*
    Float x = e[0], y = e[1], z = e[2];
    theta = acos(z);
    if (x == 0) {
        if (y >= 0) {
            phi = M_PI * 0.5;
        } else {
            phi = M_PI * 1.5;
        }
    } else {
        phi = atan2(y, x);
        if (phi < 0) phi += 2 * M_PI; // shift to [0, 2 * M_PI]
    }*/

    // y as up-vector
    Float x = e[0], y = e[1], z = e[2];
    theta = acos(y);
    if (x == 0) {
        if (z >= 0) {
            phi = A_PI * 0.5f;
        } else {
            phi = A_PI * 1.5f;
        }
    } else {
        phi = atan2(z, x);
        if (phi < 0) phi += 2 * A_PI; // shift to [0, 2 * M_PI]
    }
    
    // FIXME: should be phi += A_PI all the time as in env_light?
}

inline void Vec3::rotate_about_x(Float angle) {
    Float cosa = cos(angle);
    Float sina = sin(angle);
    Float vy =   y() * cosa + sina * z();
    Float vz = - y() * sina + cosa * z();
    set_y(vy);
    set_z(vz);
}

inline void Vec3::rotate_about_y(Float angle) {
    Float cosa = cos(angle);
    Float sina = sin(angle);
    Float vx = x() * cosa - sina * z();
    Float vz = x() * sina + cosa * z();
    set_x(vx);
    set_z(vz);
}
    
inline void Vec3::rotate_about_z(Float angle) {
    Float cosa = cos(angle);
    Float sina = sin(angle);
    Float vx = x() * cosa - sina * y();
    Float vy = x() * sina + cosa * y();
    set_x(vx);
    set_y(vy);
}

inline Vec3 rotate_about_x(const Vec3 &v, Float angle) {
    Float cosa = cos(angle);
    Float sina = sin(angle);
    return Vec3(v.x(), v.y() * cosa + sina * v.z(), - v.y() * sina + cosa * v.z());
}
    
inline Vec3 rotate_about_y(const Vec3 &v, Float angle) {
    Float cosa = cos(angle);
    Float sina = sin(angle);
    return Vec3(v.x() * cosa - sina * v.z(), v.y(), v.x() * sina + cosa * v.z());
}

inline Vec3 rotate_about_z(const Vec3 &v, Float angle) {
    Float cosa = cos(angle);
    Float sina = sin(angle);    
    return Vec3(v.x() * cosa - sina * v.y(), v.x() * sina + cosa * v.y(), v.z());
}

inline void Vec3::clamp(const Vec3 &a, const Vec3 &b) {
    e[0] = std::min(std::max(e[0], a.x()), b.x());
    e[1] = std::min(std::max(e[1], a.y()), b.y());
    e[2] = std::min(std::max(e[2], a.z()), b.z());
}

inline Vec3 clamp(const Vec3 &u, const Vec3 &a, const Vec3 &b) {
    return Vec3(
        std::min(std::max(u.x(), a.x()), b.x()),
        std::min(std::max(u.y(), a.y()), b.y()),
        std::min(std::max(u.z(), a.z()), b.z()));
}

} // end namespace
#endif

