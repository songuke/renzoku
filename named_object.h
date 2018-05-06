#ifndef _NAMED_OBJECT_H_
#define _NAMED_OBJECT_H_

#include <string>
using namespace std;

namespace Renzoku {

class NamedObject {
public:
    void set_name(string name) {
        this->name = name;
    }

    string get_name() const {
        return name;
    }

protected:
    string name;
};

} // end namespace

#endif