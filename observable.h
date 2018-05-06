#ifndef _OBSERVABLE_H_
#define _OBSERVABLE_H_

#include <vector>
using namespace std;

namespace Renzoku {

template <typename ObserverType>
class Observable {
public:
    void attach_observer(ObserverType *ob) {
        for (int i = 0; i < observers.size(); ++i) {
            if (observers[i] == ob) {
                return;
            }
        }
        observers.push_back(ob);
    }

    void detach_observer(ObserverType *ob) {
        /*
        // this code does not compile on OS X
        for (ObserverTypes::iterator i = observers.begin(); i != observers.end(); ++i) {
            if ((*i) == ob) {
                observers.erase(i);
                break;
            }
        }*/
        for (int i = 0; i < observers.size(); ++i) {
            if (observers[i] == ob) {
                observers.erase(observers.begin() + i);
                break;
            }
        }   
    }

    void clear() {
        observers.clear();
    }

    int count_observers() {
        return observers.size();
    }

protected:
    vector<ObserverType*> observers;
};

} // end namespace
#endif