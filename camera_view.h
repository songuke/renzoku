#ifndef _CAMERA_VIEW_H_
#define _CAMERA_VIEW_H_

#include "common.h"
#include "named_object.h"

namespace Renzoku {

class CameraView : public NamedObject {
public:
    CameraView() {}

    Camera *get_camera() const {
        return this->camera;
    }

    void set_camera(Camera *camera) {
        this->camera = camera;
    }

    Integrator *get_integrator() const {
        return integrator;
    }

    void set_integrator(Integrator *integrator) {
        this->integrator = integrator;
    }

private:
    Camera *camera;
    Integrator *integrator;
}; 

} // end namespace

#endif