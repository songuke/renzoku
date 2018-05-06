#ifndef _MIRROR_H_
#define _MIRROR_H_

#include "common.h"
#include "bsdf.h"
#include "rgb.h"

namespace Renzoku {

class Mirror : public Bsdf {
public:
    Mirror();
    Mirror(const Rgb &kd);
    
    Rgb sample(Random &rd, const Onb &uvn, const Vec3 &wo, Vec3 &wi, Float &pdf);    
    
    virtual Rgb eval(const Onb &uvn, const Vec3 &wo, const Vec3 &wi);  
    virtual Float pdf(const Onb &uvn, const Vec3 &wo, const Vec3 &wi);
    
    virtual bool    is_singular() const              { return true; }
    virtual Rgb     get_representative_color() const { return kd; }
    virtual Bsdf::Type get_bsdf_type() const { return Bsdf::MIRROR; }

protected:
    Rgb kd;
};

} // end namespace

#endif
