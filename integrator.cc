#include "integrator.h"

namespace Renzoku {

Integrator::Integrator() : scene(NULL), max_bounce(4), suffix("unknown") {
}

void Integrator::initialize(Scene *scene) {
    this->scene = scene;
}

inline bool Integrator::is_viewer_outputing() const {
    return false;
}

} // end namespace