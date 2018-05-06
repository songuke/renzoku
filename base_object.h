#ifndef _BASE_OBJECT_
#define _BASE_OBJECT_

#include <cstddef>
using namespace std;

namespace Renzoku {
class BaseObject
{
public:
    BaseObject(void);
    virtual ~BaseObject(void);

    void *operator new(size_t size); // ensure 16-byte alignment when new/delete is called
    void operator delete(void *ptr);
    void *operator new[](size_t size); 
    void operator delete[](void *ptr);
};

}; // end namespace Renzoku

#endif
