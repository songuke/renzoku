#ifndef _PARAM_DICT_
#define _PARAM_DICT_

#include <boost/any.hpp>

#include <string>
#include <map>
using namespace std;
using namespace boost;

#include "vec2.h"
#include "vec3.h"

namespace Renzoku {

class ParamDict {
public:
    float   get_float(const string &s);
    int     get_int(const string &s);
    bool    get_bool(const string &s);
    double  get_double(const string &s);
    string  get_string(const string &s);
    Vec3    get_vec3(const string &s);
    Vec2    get_vec2(const string &s);
    
    float   get_float(const string &s,  float default_value);
    int     get_int(const string &s,    int default_value);
    bool    get_bool(const string &s,   bool default_value);
    double  get_double(const string &s, double default_value);
    string  get_string(const string &s, const string &default_value);
    Vec3    get_vec3(const string &s,   const Vec3 &default_value);
    Vec2    get_vec2(const string &s,   const Vec2 &default_value);
        
    void    add(const string &s, float v);
    void    add(const string &s, int v);
    void    add(const string &s, bool v);
    void    add(const string &s, double v);
    void    add(const string &s, const string &v);
    void    add(const string &s, const Vec2 &v);
    void    add(const string &s, const Vec3 &v);

    void    clear();

protected:
    /**
     * Helper function. It can add any objects to the dictionary, and
     * can cause type incompatibility (add an Rgb and retrieve with get_vec3).
     * So we prefer just using it internally.
     */ 
    void    add(const string &s, const boost::any &a);

protected:
    std::map<std::string, boost::any> m;
};

} // end namespace

#endif


