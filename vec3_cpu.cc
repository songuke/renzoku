#include "vec3_cpu.h"

namespace Renzoku {

istream& operator>>(istream &is, Vec3 &u) {
    is >> u.e[0]; // x
    is >> u.e[1]; // y
    is >> u.e[2]; // z    
    return is;
}

ostream& operator<<(ostream &os, const Vec3 &u) {         
    os << u.e[0] << ", " << u.e[1] << ", " << u.e[2]; 
    return os;
}
}; // end namespace Rzrt
