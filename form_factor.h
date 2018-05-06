#ifndef _FORM_FACTOR_H_
#define _FORM_FACTOR_H_

#include "vec3.h"

namespace Renzoku {

class FormFactor {

public:
    static inline Float form_factor_abs(const Vec3 &p, const Vec3 &pn, 
                                        const Vec3 &q, const Vec3 &qn) {

        Vec3 pq = q - p;
        Float dot1 =   dot(pn, pq);
        Float dot2 = - dot(qn, pq);    
        Float dot_pq = dot(pq, pq);
        Float form_factor = (dot1 * dot2) / (dot_pq * dot_pq);    
        return fabs(form_factor);
    }

    static inline Float form_factor(const Vec3 &p, const Vec3 &pn, 
                                    const Vec3 &q, const Vec3 &qn) {

        Vec3 pq = q - p;
        Float dot1 =   dot(pn, pq);
        Float dot2 = - dot(qn, pq);    
        Float dot_pq = dot(pq, pq);
        Float form_factor = (dot1 * dot2) / (dot_pq * dot_pq);
        return std::max(0.0f, form_factor);
    }
        
};

}

#endif