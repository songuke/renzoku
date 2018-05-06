#ifndef _LINEAR_H_
#define _LINEAR_H_

#include "common.h"

namespace Renzoku {

class Linear {
public:
    static inline void 
        solve33(Float a, Float b, Float c,
                Float d, Float e, Float f,
                Float g, Float h, Float i,
                Float p, Float q, Float r,
                Float &x, Float &y, Float &z);

    static inline Float 
        det(Float a, Float b, Float c,
            Float d, Float e, Float f,
            Float g, Float h, Float i);
};

inline void Linear::solve33(Float a, Float b, Float c,
                            Float d, Float e, Float f,
                            Float g, Float h, Float i,
                            Float p, Float q, Float r,
                            Float &x, Float &y, Float &z) {

    Float inv_denom = 1.0f / det(a, b, c, d, e, f, g, h , i);
    x = inv_denom * det(p, b, c, 
                        q, e, f, 
                        r, h, i);

    y = inv_denom * det(a, p, c, 
                        d, q, f, 
                        g, r, i);

    z = inv_denom * det(a, b, p, 
                        d, e, q, 
                        g, h, r);
}

inline Float Linear::det(Float a, Float b, Float c,
                  Float d, Float e, Float f,
                  Float g, Float h, Float i) {
    return (a * e * i + b * f * g + c * d * h) - (c * e * g + b * d * i + a * f * h);
}

} // end namespace

#endif