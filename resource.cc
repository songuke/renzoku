#include "resource.h"

namespace Renzoku {

ID Resource::next_id = 0;

ID Resource::get_next_id() {
    ID next = next_id;
    ++next_id;
    return next;
}

} // end namespace Renzoku