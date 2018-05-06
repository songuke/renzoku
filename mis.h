#ifndef _MIS_H_
#define _MIS_H_

#include "common.h"

namespace Renzoku {

class Mis {

public:
    static inline Float balance(Float a, Float b) {
        return a / (a + b);
    }

    static inline Float power(Float a, Float b) {
        return (a * a) / (a * a + b * b);
    }
};

} 

#endif