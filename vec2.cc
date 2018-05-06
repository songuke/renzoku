#include "vec2.h"

namespace Renzoku {
istream& operator>>(istream &is, Vec2 &v) {
    Float tmp;
    is >> tmp; v.e[0] = tmp;
    is >> tmp; v.e[1] = tmp;
    return is;
}

ostream& operator<<(ostream &os, const Vec2 &v) {
    os << v.e[0] << ' ' << v.e[1];
    return os;
}
} // end namespace Renzoku
