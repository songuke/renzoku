#ifndef _RESOURCE_H_
#define _RESOURCE_H_

#include "common.h"

namespace Renzoku {

/** 
 * Resource management.
 */
class Resource {
public:
    static ID get_next_id();

protected:
    static ID next_id;
};

} // end namespace Renzoku

#endif
