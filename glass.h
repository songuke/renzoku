#ifndef _GLASS_H_
#define _GLASS_H_

#include "common.h"
#include "bsdf.h"
#include "rgb.h"

namespace Renzoku {

class Glass : public Bsdf {
public:
    Glass(Float inner);
    Glass(Float inner, Float outer);
    Glass(Float inner, const Rgb &kd);
    Glass(Float inner, Float outer, const Rgb &kd);
    

    virtual Rgb sample(Random &rd, const Onb &uvn, const Vec3 &wo, Vec3 &wi, Float &pdf);
    virtual Rgb eval(const Onb &uvn, const Vec3 &wo, const Vec3 &wi);
    virtual Float pdf(const Onb &uvn, const Vec3 &wo, const Vec3 &wi);

    virtual bool is_singular()     const { return true; }
    virtual bool is_transmissive() const { return true; }
    Rgb get_representative_color() const { return kd; }

    virtual Bsdf::Type get_bsdf_type() const { return Bsdf::GLASS; }

    inline Float get_inner_index() const;
    inline Float get_outer_index() const;

private:
    Rgb sample(Random &rd, const Vec3 &nn, Float nr, Float nt, const Vec3 &wo, Vec3 &wi, Float &pdf);

protected:

    Float inner, outer;     // medium index. Outer aligns with the surface normal, and inner aligns with reverse normal.
    Rgb kd;
};

inline Float Glass::get_inner_index() const {
    return inner;
}

inline Float Glass::get_outer_index() const {
    return outer;
}

} // end namespace
#endif

