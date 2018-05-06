#include "material.h"
#include "lambertian.h"

namespace Renzoku {

BsdfMaterial *DefaultMaterial::_white = new BsdfMaterial(new Lambertian(1.0f, 1.0f, 1.0f));
BsdfMaterial *DefaultMaterial::_pi    = new BsdfMaterial(new Lambertian(A_PI, A_PI, A_PI));

BsdfMaterial *DefaultMaterial::white() {
	return _white;
}

BsdfMaterial *DefaultMaterial::pi() {
    return _pi;
}

} // end namespace
