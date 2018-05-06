#ifndef _NUMBER_H_
#define _NUMBER_H_

#include "common.h"

namespace Renzoku {

inline bool is_nan(Float val) {
    return (val != val) ? true : false;
}

} // end namespace

#endif