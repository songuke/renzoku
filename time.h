#ifndef _TIME_H_
#define _TIME_H_

#include <string>
using namespace std;

#include "boost/date_time/posix_time/posix_time.hpp"
using namespace boost::posix_time;

namespace Renzoku {

class Time {
public:
    static string get_system_time_string() {
        ptime now(microsec_clock::local_time());
        return to_simple_string(now);
    }
};

} // end namespace

#endif