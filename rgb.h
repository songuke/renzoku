#ifndef _RGB_H_
#define _RGB_H_

#include <cmath>
#include <iostream>
using namespace std;

#include "common.h"
#include "vec3.h"

namespace Renzoku {


class Rgb {
public: 
    Rgb() : r(0), g(0), b(0) {}
    Rgb(Float t) : r(t), g(t), b(t) {}
    Rgb(Float _r, Float _g, Float _b) : r(_r), g(_g), b(_b) {}
    Rgb(const Rgb &c) : r(c.r), g(c.g), b(c.b) {}
    Rgb(Vec3 v) : r(v.x()), g(v.y()), b(v.z()) {}

    void set_red(Float _r) { r = _r; }
    void set_green(Float _g) { g = _g; }
    void set_blue(Float _b) { b = _b; }
    
    bool is_nan() const { return (r != r || g != g || b != b); }

    Rgb& operator=(const Rgb& c);
    Rgb& operator=(Float t);
    Rgb& operator+=(const Rgb& c);
    Rgb& operator*=(const Rgb& c);
    Rgb& operator/=(const Rgb& c);
    Rgb& operator*=(Float t);
    Rgb& operator/=(Float t);

    Rgb& operator+() { return *this; }
    Rgb operator-() const { return Rgb(-r, -g, -b); }
    
    Float red() const { return r; }
    Float green() const { return g; }
    Float blue() const { return b; }
    
    Float avg() const { return (r + g + b) / 3.0f; }
    Float max() const { 
        Float tmp = r;
        if (g > tmp) tmp = g;
        if (b > tmp) tmp = b;
        return tmp;
    }
    Float min() const {
        Float tmp = r;
        if (g < tmp) tmp = g;
        if (b < tmp) tmp = b;
        return tmp;
    }

    Float luminance() const { return sqrt(r*r + g*g + b*b); }
    
    void clamp();
    void clamp(Float threshold);

    /**
     * Convert from HSV color space to RGB.
     */
    static Rgb from_hsv(Float hue, Float saturation, Float value); 
    
    friend ostream& operator<<(ostream &os, const Rgb& c);
    friend Rgb operator+(const Rgb& c1, const Rgb& c2);
    friend Rgb operator-(const Rgb& c1, const Rgb& c2);
    friend Rgb operator*(const Rgb& c1, const Rgb& c2);
    friend Rgb operator/(const Rgb& c1, const Rgb& c2);
    friend Rgb operator/(const Rgb& c, Float t);
    friend Rgb operator*(const Rgb& c, Float t);
    friend Rgb operator*(Float t, const Rgb& c);
 
    friend bool operator==(const Rgb &v1, const Rgb &v2);
    friend bool operator!=(const Rgb &v1, const Rgb &v2);   
    friend bool operator<(const Rgb &v1, const Rgb &v2);
    friend bool operator<=(const Rgb &v1, const Rgb &v2);

    /**
     * V component in HSV.
     */
    friend Float value(const Rgb &c);    
    Float value() const { return sqrt(r*r + g*g + b*b); }
    
    friend Rgb abs(const Rgb &c);

    Vec3 to_vec3() const;

public:
    Float r, g, b; 
};

typedef vector<Rgb> Rgbs;

class DefaultRgb {
public:
    static const Rgb white;
    static const Rgb black;
    static const Rgb red;
    static const Rgb green;
    static const Rgb blue;
    static const Rgb yellow;
    static const Rgb cyan;
    static const Rgb magenta;
    static const Rgb grey;
};

inline Rgb& Rgb::operator=(const Rgb &c) {
    r = c.r; 
    g = c.g;
    b = c.b;
    return *this;
}

inline Rgb& Rgb::operator=(Float t) {
    r = g = b = t;
    return *this;
}

inline Rgb& Rgb::operator+=(const Rgb& c) {
    return *this = *this + c; }

inline Rgb& Rgb::operator*=(const Rgb& c) {
    return *this = *this * c; }

inline Rgb& Rgb::operator/=(const Rgb& c) {
    return *this = *this / c; }

inline Rgb& Rgb::operator*=(Float c) {
    return *this = *this * c; }

inline Rgb& Rgb::operator/=(Float c) {
    return *this = *this / c; }

inline void Rgb::clamp() {    
    r = r < 0 ? 0 : (r > 1 ? 1 : r);
    g = g < 0 ? 0 : (g > 1 ? 1 : g);
    b = b < 0 ? 0 : (b > 1 ? 1 : b);
}

inline void Rgb::clamp(Float threshold) {
    r = r > threshold ? threshold : r;
    g = g > threshold ? threshold : g;
    b = b > threshold ? threshold : b;
}

inline ostream& operator<<(ostream &os, const Rgb& c) {
    os << c.r << ' ' << c.g << ' ' << c.b;
    return os;
}

inline Rgb operator+(const Rgb& c1, const Rgb& c2) {
    return Rgb(c1.r + c2.r, c1.g + c2.g, c1.b + c2.b); }

inline Rgb operator-(const Rgb& c1, const Rgb& c2) {
    return Rgb(c1.r - c2.r, c1.g - c2.g, c1.b - c2.b); }

inline Rgb operator*(const Rgb& c1, const Rgb& c2) {
    return Rgb(c1.r * c2.r, c1.g * c2.g, c1.b * c2.b); }

inline Rgb operator/(const Rgb& c1, const Rgb& c2) {
    return Rgb(c1.r / c2.r, c1.g / c2.g, c1.b / c2.b); }

inline Rgb operator/(const Rgb& c, Float t) {
    return Rgb(c.r / t, c.g / t, c.b / t); }

inline Rgb operator*(const Rgb& c, Float t) {
    return Rgb(c.r * t, c.g * t, c.b * t); }

inline Rgb operator*(Float t, const Rgb& c) {
    return Rgb(c.r * t, c.g * t, c.b * t); }

inline bool operator==(const Rgb &v1, const Rgb &v2) {
    return ((v1.r == v2.r) && (v1.g == v2.g) && (v1.b == v2.b));
}

inline bool operator!=(const Rgb &v1, const Rgb &v2) {
    return ((v1.r != v2.r) || (v1.g != v2.g) || (v1.b != v2.b));
}

inline bool operator<(const Rgb &v1, const Rgb &v2) {
    return ((v1.r < v2.r) && (v1.g < v2.g) && (v1.b < v2.b));
}

inline bool operator<=(const Rgb &v1, const Rgb &v2) {
    return ((v1.r <= v2.r) && (v1.g <= v2.g) && (v1.b <= v2.b));
}

inline Vec3 Rgb::to_vec3() const {
    return Vec3(r, g, b); 
}

inline Float value(const Rgb &c) {
    return sqrt(c.r * c.r + c.g * c.g + c.b * c.b);
}

inline Rgb abs(const Rgb &c) {
    return Rgb(fabs(c.r), fabs(c.g), fabs(c.b));
}

/*
class Rgb {
public: 
    Rgb() { e[0] = 0; e[1] = 0; e[2] = 0; }
    Rgb(Float _r, Float _g, Float _b) { e[0] = _r; e[1] = _g; e[2] = _b;}
    Rgb(const Rgb &c) { e[0] = c.e[0]; e[1] = c.e[1]; e[2] = c.e[2]; }

    void set_red(Float _r) { e[0] = _r; }
    void set_green(Float _g) { e[1] = _g; }
    void set_blue(Float _b) { e[2] = _b; }
    
    Rgb& operator=(const Rgb& c);
    Rgb& operator+=(const Rgb& c);
    Rgb& operator*=(const Rgb& c);
    Rgb& operator/=(const Rgb& c);
    Rgb& operator*=(Float t);
    Rgb& operator/=(Float t);

    Rgb& operator+() { return *this; }
    Rgb operator-() const { return Rgb(-e[0], -e[1], -e[2]); }
    
    Float red() const { return e[0]; }
    Float green() const { return e[1]; }
    Float blue() const { return e[2]; }
    
    Float avg() const { return (e[0] + e[1] + e[2]) / 3.0; }
    Float max() const { 
        Float tmp = e[0];
        if (e[1] > tmp) tmp = e[1];
        if (e[2] > tmp) tmp = e[2];
        return tmp;
    }
    Float min() const {
        Float tmp = e[0];
        if (e[1] < tmp) tmp = e[1];
        if (e[2] < tmp) tmp = e[2];
        return tmp;
    }

    Float luminance() const { return sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]); }
    void clamp();
 
    friend ostream& operator<<(ostream &os, const Rgb& c);
    friend Rgb operator+(const Rgb& c1, const Rgb& c2);
    friend Rgb operator*(const Rgb& c1, const Rgb& c2);
    friend Rgb operator/(const Rgb& c1, const Rgb& c2);
    friend Rgb operator/(const Rgb& c, Float t);
    friend Rgb operator*(const Rgb& c, Float t);
    friend Rgb operator*(Float t, const Rgb& c);
 
    friend bool operator==(const Rgb &v1, const Rgb &v2);
    friend bool operator!=(const Rgb &v1, const Rgb &v2);   
public:
    Float e[3];
};

inline Rgb& Rgb::operator=(const Rgb &c) {
    e[0] = c.e[0];
    e[1] = c.e[1];
    e[2] = c.e[2];
    return *this;
}

inline Rgb& Rgb::operator+=(const Rgb& c) {
    return *this = *this + c; }

inline Rgb& Rgb::operator*=(const Rgb& c) {
    return *this = *this * c; }

inline Rgb& Rgb::operator/=(const Rgb& c) {
    return *this = *this / c; }

inline Rgb& Rgb::operator*=(Float c) {
    return *this = *this * c; }

inline Rgb& Rgb::operator/=(Float c) {
    return *this = *this / c; }

inline void Rgb::clamp() {
    if (e[0] < 0.) e[0] = 0.;
    if (e[0] > 1.) e[0] = 1.;
    if (e[1] < 0.) e[1] = 0.;
    if (e[1] > 1.) e[1] = 1.;
    if (e[2] < 0.) e[2] = 0.;
    if (e[2] > 1.) e[2] = 1.;
}

inline ostream& operator<<(ostream &os, const Rgb& c) {
    os << c.e[0] << ' ' << c.e[1] << ' ' << c.e[2];
    return os;
}

inline Rgb operator+(const Rgb& c1, const Rgb& c2) {
    return Rgb(c1.e[0] + c2.e[0], c1.e[1] + c2.e[1], c1.e[2] + c2.e[2]); }

inline Rgb operator*(const Rgb& c1, const Rgb& c2) {
    return Rgb(c1.e[0] * c2.e[0], c1.e[1] * c2.e[1], c1.e[2] * c2.e[2]); }

inline Rgb operator/(const Rgb& c1, const Rgb& c2) {
    return Rgb(c1.e[0] / c2.e[0], c1.e[1] / c2.e[1], c1.e[2] / c2.e[2]); }

inline Rgb operator/(const Rgb& c, Float t) {
    return Rgb(c.e[0] / t, c.e[1] / t, c.e[2] / t); }

inline Rgb operator*(const Rgb& c, Float t) {
    return Rgb(c.e[0] * t, c.e[1] * t, c.e[2] * t); }

inline Rgb operator*(Float t, const Rgb& c) {
    return Rgb(c.e[0] * t, c.e[1] * t, c.e[2] * t); }

inline bool operator==(const Rgb &v1, const Rgb &v2) {
    if (v1.e[0] != v2.e[0]) return false;
    if (v1.e[1] != v2.e[1]) return false;
    if (v1.e[2] != v2.e[2]) return false;
    return true;
}

inline bool operator!=(const Rgb &v1, const Rgb &v2) {
    return !(v1 == v2);
}
*/
}; // end namespace Renzoku
#endif
