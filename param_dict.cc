#include "param_dict.h"

namespace Renzoku {

float   ParamDict::get_float(const string &s) {
    try {
        return boost::any_cast<float>(m[s]);
    } catch (std::exception e) {        
        return (float)boost::any_cast<int>(m[s]);
    }
}

int     ParamDict::get_int(const string &s) {
    // TODO: check if s is available before m[s] and throw more informative exception
    return boost::any_cast<int>(m[s]);
}

bool    ParamDict::get_bool(const string &s) {
    // TODO: check if s is available before m[s] and throw more informative exception
    return boost::any_cast<int>(m[s]);
}

double  ParamDict::get_double(const string &s) {
    return boost::any_cast<double>(m[s]);
}

string  ParamDict::get_string(const string &s) {
    return boost::any_cast<string>(m[s]);
}

Vec3    ParamDict::get_vec3(const string &s) {
    return boost::any_cast<Vec3>(m[s]);
}

Vec2    ParamDict::get_vec2(const string &s) {
    return boost::any_cast<Vec2>(m[s]);
}

float   ParamDict::get_float(const string &s,  float default_value) {    
    return (m.find(s) == m.end()) ? default_value : get_float(s);
}

int     ParamDict::get_int(const string &s,    int default_value) {
    return (m.find(s) == m.end()) ? default_value : get_int(s);
}

bool    ParamDict::get_bool(const string &s, bool default_value) {
    return (m.find(s) == m.end()) ? default_value : get_bool(s);
}

double  ParamDict::get_double(const string &s, double default_value) {
    return (m.find(s) == m.end()) ? default_value : get_double(s);
}

string  ParamDict::get_string(const string &s, const string &default_value) {
    return (m.find(s) == m.end()) ? default_value : get_string(s);
}

Vec3    ParamDict::get_vec3(const string &s,   const Vec3 &default_value) {
    return (m.find(s) == m.end()) ? default_value : get_vec3(s);
}

Vec2    ParamDict::get_vec2(const string &s,   const Vec2 &default_value) {
    return (m.find(s) == m.end()) ? default_value : get_vec2(s);
}

void    ParamDict::add(const string &s, const boost::any &a) {
    m[s] = a;
}
    
void    ParamDict::add(const string &s, float v) {    
    add(s, boost::any(v));
}

void    ParamDict::add(const string &s, int v) {
    add(s, boost::any(v));
}

void    ParamDict::add(const string &s, bool v) {
    add(s, boost::any(v));
}

void    ParamDict::add(const string &s, double v) {
    add(s, boost::any(v));
}

void    ParamDict::add(const string &s, const string &v) {
    add(s, boost::any(v));
}

void    ParamDict::add(const string &s, const Vec2 &v) {
    add(s, boost::any(v));
}

void    ParamDict::add(const string &s, const Vec3 &v) {
    add(s, boost::any(v));
}

void    ParamDict::clear() {
    m.clear();
}

} // end namespace