#include "thread_controller.h"

namespace Renzoku {

void ThreadController::add_thread(boost::thread *thread) {
    threads.push_back(thread);
}

void ThreadController::kill_all() {
    for (int i = 0; i < threads.size(); ++i) {
        threads[i]->interrupt();
        threads[i]->join();
    }
    threads.clear();
}

} // end namespace