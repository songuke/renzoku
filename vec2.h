#ifndef _VEC2_H_
#define _VEC2_H_

#include <cmath>
#include <iostream>
using namespace std;

#include "common.h"
#include "random.h"

namespace Renzoku {
class Vec2 {
public:
    Vec2();
    Vec2(Float v);
    Vec2(Float e0, Float e1);
    Vec2(const Vec2 &v) { *this = v; }
    
    Float x() const { return e[0]; }
    Float y() const { return e[1]; }
    
    const Vec2& operator+() const;
    Vec2 operator-() const; 
    
    Float operator[](int i) const { return e[i]; }
    Float& operator[](int i) { return e[i]; }

    Float length() const;
    Float squared_length() const;
    
    void make_unit_vector();

    void set_x(Float _x) { e[0] = _x; }
    void set_y(Float _y) { e[1] = _y; }

    Float min_component() const;
    Float max_component() const;
    Float min_abs_component() const;
    Float max_abs_component() const;

    int arg_min_component() const;
    int arg_max_component() const;
    int arg_min_abs_component() const;
    int arg_max_abs_component() const;

    friend bool operator==(const Vec2 &v1, const Vec2 &v2);
    friend bool operator!=(const Vec2 &v1, const Vec2 &v2);
    
    friend istream& operator>>(istream &is, Vec2 &t);
    friend ostream& operator<<(ostream &os, const Vec2 &t);
  
    friend Vec2 operator+(const Vec2 &v1, const Vec2 &v2);
    friend Vec2 operator-(const Vec2 &v1, const Vec2 &v2);
    friend Vec2 operator*(const Vec2 &v1, const Vec2 &v2);
    friend Vec2 operator/(const Vec2 &v1, const Vec2 &v2);

    friend Vec2 operator-(const Vec2 &v, Float t);
    friend Vec2 operator-(Float t, const Vec2 &v);
    friend Vec2 operator/(const Vec2 &v, Float t);
    friend Vec2 operator*(const Vec2 &v, Float t);
    friend Vec2 operator*(Float scalar, const Vec2 &v);
    
    
    Vec2& operator=(const Vec2 &v);
    Vec2& operator+=(const Vec2 &v);
    Vec2& operator-=(const Vec2 &v);
    Vec2& operator*=(const Float t);
    Vec2& operator/=(const Float t);

    friend Vec2 unit_vector(const Vec2 &v);
    friend Vec2 min(const Vec2 &v1, const Vec2 &v2);
    friend Vec2 max(const Vec2 &v1, const Vec2 &v2);
    friend Float dot(const Vec2 &v1, const Vec2 &v2);

    void random(Random &rd) { e[0] = rd(); e[1] = rd(); }

public:
    Float e[2];
};

inline Vec2::Vec2() {
    e[0] = 0;
    e[1] = 0;
}

inline Vec2::Vec2(Float v) {
    e[0] = e[1] = v;
}

inline Vec2::Vec2(Float e0, Float e1) { 
    e[0] = e0; e[1] = e1; }

inline const Vec2& Vec2::operator+() const { return *this; }

inline Vec2 Vec2::operator-() const { return Vec2(-e[0], -e[1]); }

inline Float Vec2::length() const { 
    return sqrt(e[0]*e[0] + e[1]*e[1]); }

inline Float Vec2::squared_length() const {
    return e[0]*e[0] + e[1]*e[1]; }

inline void Vec2::make_unit_vector() {
    *this = *this / (*this).length(); }

inline Float Vec2::min_component() const {
    Float tmp = e[0];
    if (e[1] < tmp) tmp = e[1];
    return tmp; 
}

inline Float Vec2::max_component() const {
    Float tmp = e[0];
    if (e[1] > tmp) tmp = e[1];
    return tmp;
}

inline Float Vec2::min_abs_component() const {
    Float tmp = fabs(e[0]);
    if (fabs(e[1]) < tmp) tmp = fabs(e[1]);
    return tmp;
}

inline Float Vec2::max_abs_component() const {
    Float tmp = fabs(e[0]);
    if (fabs(e[1]) > tmp) tmp = fabs(e[1]);
    return tmp;
}

inline int Vec2::arg_min_component() const {
    return (e[0] < e[1]) ? 0 : 1;
}

inline int Vec2::arg_max_component() const {
    return (e[0] > e[1]) ? 0 : 1;
}

inline int Vec2::arg_min_abs_component() const {
    return (fabs(e[0]) < fabs(e[1])) ? 0 : 1;
}

inline int Vec2::arg_max_abs_component() const {
    return (fabs(e[0]) > fabs(e[1])) ? 0 : 1;
}

inline bool operator==(const Vec2 &v1, const Vec2 &v2) {
    if (v1.e[0] != v2.e[0]) return false;
    if (v1.e[1] != v2.e[1]) return false;
    return true;
}

inline bool operator!=(const Vec2 &v1, const Vec2 &v2) {
    return !(v1 == v2); }

inline Vec2 operator+(const Vec2 &v1, const Vec2 &v2) {
    return Vec2(v1.e[0] + v2.e[0], v1.e[1] + v2.e[1]); 
}

inline Vec2 operator-(const Vec2 &v1, const Vec2 &v2) {
    return Vec2(v1.e[0] - v2.e[0], v1.e[1] - v2.e[1]);
}

inline Vec2 operator*(const Vec2 &v1, const Vec2 &v2) {
    return Vec2(v1.e[0] * v2.e[0], v1.e[1] * v2.e[1]); 
}

inline Vec2 operator/(const Vec2 &v1, const Vec2 &v2) {
    return Vec2(v1.e[0] / v2.e[0], v1.e[1] / v2.e[1]); 
}

inline Vec2 operator-(const Vec2 &v, Float t) {
    return Vec2(v.e[0] - t, v.e[1] - t); }

inline Vec2 operator-(Float t, const Vec2 &v) {
    return Vec2(v.e[0] - t, v.e[1] - t); } 

inline Vec2 operator/(const Vec2 &v, Float t) {
    return Vec2(v.e[0]/t, v.e[1]/t); }

inline Vec2 operator*(const Vec2 &v, Float t) {
    return Vec2(v.e[0]*t, v.e[1]*t); }

inline Vec2 operator*(Float t, const Vec2 &v) {
    return Vec2(v.e[0]*t, v.e[1]*t); }

inline Vec2& Vec2::operator=(const Vec2 &v) {
    e[0] = v.e[0]; e[1] = v.e[1]; return *this; }

inline Vec2& Vec2::operator+=(const Vec2 &v) {
    *this = *this + v; return *this; }

inline Vec2& Vec2::operator-=(const Vec2 &v) {
    *this = *this - v; return *this; }

inline Vec2& Vec2::operator*=(Float t) {
    *this = *this * t; return *this; }

inline Vec2& Vec2::operator/=(Float t) {
    *this = *this / t; return *this; }

inline Float dot(const Vec2 &v1, const Vec2 &v2) {
    return v1.e[0]*v2.e[0] + v1.e[1]*v2.e[1]; }

inline Vec2 unit_vector(const Vec2 &v) {
    return v / v.length(); }

inline Vec2 min(const Vec2 &v1, const Vec2 &v2) {
    Vec2 tmp(v1);
    if (v2.x() < v1.x()) tmp.set_x(v2.x());
    if (v2.y() < v1.y()) tmp.set_y(v2.y());
    return tmp;
}

inline Vec2 max(const Vec2 &v1, const Vec2 &v2) {
    Vec2 tmp(v1);
    if (v2.x() > v1.x()) tmp.set_x(v2.x());
    if (v2.y() > v1.y()) tmp.set_y(v2.y());
    return tmp;
}

} // end namespace Renzoku
#endif

