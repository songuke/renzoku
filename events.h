#ifndef _EVENTS_H_
#define _EVENTS_H_

// this file is shared between trackball C code and other C++ files. 

#ifdef __cplusplus
    extern "C" {
#endif

typedef enum {
    MOUSE_BUTTON_LEFT,
    MOUSE_BUTTON_RIGHT,
    MOUSE_BUTTON_MIDDLE
} MouseButton;

typedef enum {
    MOUSE_STATE_DOWN,
    MOUSE_STATE_UP,
    MOUSE_STATE_DOUBLE
} MouseState;

#ifdef __cplusplus
    }
#endif

#ifdef __cplusplus

#include "observable.h"

#include <vector>
using namespace std;

namespace Renzoku {

class IMouseObserver {
public:
    virtual void on_click(MouseButton button, MouseState state, int x, int y) = 0;
};

class MouseObservable : public Observable<IMouseObserver> {
public:

    void notify_click(MouseButton button, MouseState state, int x, int y) {
        for (int i = 0; i < observers.size(); ++i)
            observers[i]->on_click(button, state, x, y);
    }

};

}   // end namespace

#endif


#endif
