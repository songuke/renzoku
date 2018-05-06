#ifndef _THREAD_CONTROLLER_H_
#define _THREAD_CONTROLLER_H_

#include "common.h"
#include <boost/thread.hpp>

#include <vector>
using namespace std;

namespace Renzoku {

class ThreadController {
public:
    void add_thread(boost::thread *thread);
    void kill_all();

private:
    vector<boost::thread *> threads;
};

} // end namespace

#endif