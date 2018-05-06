#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

#include "common.h"
#include "math3.h"
#include "observable.h"

namespace Renzoku {

class IntegratorObserver {
public:
    virtual void on_update_border(int tid, int x0, int y0, int x1, int y1) {}
    virtual void on_update_image(int tid, ImageFloat *img, int x0, int y0, int x1, int y1) = 0;
    virtual void on_complete_frame(int tid, ImageFloat *img) {}
};
typedef vector<IntegratorObserver *> IntegratorObservers;


class IntegratorObservable : public Observable<IntegratorObserver> {
public:
    void notify_update_image(int tid, ImageFloat *img, int x0, int y0, int x1, int y1) {
        IntegratorObservers::iterator i;
        for (i = observers.begin(); i != observers.end(); ++i) 
            (*i)->on_update_image(tid, img, x0, y0, x1, y1);
    }

    void notify_complete_frame(int tid, ImageFloat *img) {
        IntegratorObservers::iterator i;
        for (i = observers.begin(); i != observers.end(); ++i) 
            (*i)->on_complete_frame(tid, img);
    }
    
    void notify_update_border(int tid, int x0, int y0, int x1, int y1) {
        IntegratorObservers::iterator i;
        for (i = observers.begin(); i != observers.end(); ++i) 
            (*i)->on_update_border(tid, x0, y0, x1, y1);
    }
};

class IDirectIndirect {
public:
    IDirectIndirect() : is_direct(true), is_indirect(true) {
    }

    inline void set_direct_lighting(bool direct)         { this->is_direct = direct; }
    inline void set_indirect_lighting(bool indirect)     { this->is_indirect = indirect; }
    inline bool is_direct_lighting() const               { return this->is_direct; }
    inline bool is_indirect_lighting() const             { return this->is_indirect; }

protected:
    bool is_direct, is_indirect;
};

class Integrator : public IDirectIndirect {
public:
    /**
     * Do not add any arguments to the constructor, since the integrator is inherited by many classes. 
     * This helps the subclasses to have as simple constructors as possible and no need to redefine
     * the constructor with arguments.
     */
    Integrator();

    virtual void initialize(Scene *scene); 
    virtual void integrate() {}
    virtual bool integrate_partial() { integrate(); return false; }

    inline int get_max_bounce() const;
    inline void set_max_bounce(int max_bounce);

    inline virtual string get_suffix() const;
    inline virtual bool is_viewer_outputing() const;

    inline Scene *get_scene() const;

    /**
     * Copy result into Frame. 
     * CPU methods by default do nothing. 
     * GPU methods needs to download to CPU the result from GPU memory.
     */
    virtual void transfer_to_frame() {}

protected:
    Scene *scene;
    int max_bounce;
    string suffix;

public:
    IntegratorObservable observable;
};

inline Scene *Integrator::get_scene() const {
    return scene;
}

inline int Integrator::get_max_bounce() const {
    return max_bounce;
}

inline void Integrator::set_max_bounce(int max_bounce) {
    this->max_bounce = max_bounce;
}

inline string Integrator::get_suffix() const {
    return suffix;
}

} // end namespace

#endif