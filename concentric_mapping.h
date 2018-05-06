#ifndef _CONCENTRIC_MAPPING_H_
#define _CONCENTRIC_MAPPING_H_

#include "common.h"
#include "vec2.h"
#include "vec3.h"
#include "log.h"

namespace Renzoku {

/**
 * Unit square is [0, 1]. Concentric disk and square are in [-1, 1].
 * 
 * Note that there is a scale between unit square and [-1, 1] in the mapping.
 */
static Vec2 square_to_disk(Vec2 s) {

    Float phi, r;
    Float a = 2.0f * s.x() - 1.0f;    // (a,b) is now on [-1,1]^2
    Float b = 2.0f * s.y() - 1.0f;

    if (a > -b) {               // region 1 or 2
        if (a > b) {            // region 1, also |a| > |b|
            r = a;
            phi = (A_PI/4 ) * (b/a);
        } else {                // region 2, also |b| > |a|
            r = b;
            phi = (A_PI/4) * (2 - (a/b));
        }
    } else {                    // region 3 or 4
        if (a < b) {            // region 3, also |a| >= |b|, a != 0
            r = -a;
            phi = (A_PI/4) * (4 + (b/a));
        } else {                // region 4, |b| >= |a|, but a==0 and b==0 could occur.
            r = -b;
            if (b != 0.0f)
                phi = (A_PI/4) * (6 - (a/b));
            else
                phi = 0.0f;
        }
    }

    Float u = r * cos( phi );
    Float v = r * sin( phi );
    return Vec2( u, v );
}

static Vec2 disk_to_square(Vec2 d) {

    Float r = sqrt(d.x() * d.x() + d.y() * d.y());
    Float phi = atan2(d.y(), d.x());
    if (phi < -A_PI/4) phi += TWO_PI;       // in range [-pi/4,7pi/4]

    Float a, b, x, y;

    if (phi < A_PI/4) {                     // region 1
        a = r;
        b = phi * a / (A_PI/4);
    } else if (phi < 3*A_PI/4 ) {           // region 2
        b = r;
        a = -(phi - A_PI/2) * b / (A_PI/4);
    } else if (phi < 5*A_PI/4) {            // region 3
        a = -r;
        b = (phi - A_PI) * a / (A_PI/4);
    } else {                                // region 4
        b = -r;
        a = -(phi - 3*A_PI/2) * b / (A_PI/4);
    }

    x = a * 0.5f + 0.5f;                    // map to [0, 1]
    y = b * 0.5f + 0.5f;

    return Vec2(x, y);
}

/**
 * This mapping ensures uniform points on hemisphere becomes uniform points on a disk.
 *
 * See Shirley and Chiu's paper.
 *
 * Simple projection that only takes x and y results in a non-uniform distribution on disk so we want to avoid.
 */
static Vec2 uniform_hemi_to_uniform_disk(const Vec3 &wi) {
    //return Vec2(wi.x(), wi.y());

    Float r2 = 1.0f - wi.z();
    return Vec2(wi.x() / sqrt(2.0f - r2),
                wi.y() / sqrt(2.0f - r2));    
}

static Vec3 uniform_disk_to_uniform_hemi(const Vec2 &d) {
    //return Vec3(d.x(), d.y(), sqrt(1.0f - d.squared_length()));

    Float r2 = d.x() * d.x() + d.y() * d.y();
    return Vec3(d.x() * sqrt(2.0f - r2),
                d.y() * sqrt(2.0f - r2),
                1.0f - r2);
}

static Vec2 direction_to_square(const Vec3 &wi) {    
    if (wi.z() <= 0.0f) {
        Log::info() << "Only upper hemisphere is considered." << endn;
    }

    Vec2 d = uniform_hemi_to_uniform_disk(wi);
    Vec2 s = disk_to_square(d);     // square [0, 1]^2
    return s;
}

static Vec3 square_to_direction(const Vec2 &s) {
    Vec2 d = square_to_disk(s);
    Vec3 wi = uniform_disk_to_uniform_hemi(d);
    return wi;
}

}

#endif