#ifndef _MARSHAL_H_
#define _MARSHAL_H_

#include "param_dict.h"

namespace Renzoku {

/**
 * A general object to store information of all objects from different class.
 */
class MarshalObject {
public:
    ParamDict dict;
};

class IMarshalable {
public:
    virtual void marshal(MarshalObject &o) = 0;

};

} // end namespace

#endif
