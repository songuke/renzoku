#ifndef _BOUNDINGBOX_H_
#define _BOUNDINGBOX_H_

#include "vec3.h"
#include "ray.h"
#include "onb.h"

#include "log.h"

#include <xmmintrin.h>
#include <cfloat>
using namespace std;

namespace Renzoku {
    
class BoundingBox {
public:
    BoundingBox() : v_min(Vec3(FLT_MAX)), v_max(Vec3(-FLT_MAX)) {}
	BoundingBox(const Vec3 &p) : v_min(p), v_max(p) {}

    BoundingBox(const Vec3 &p0, const Vec3 &p1) {
        v_min = min(p0, p1);
        v_max = max(p0, p1);
    }
    
	BoundingBox(const Vec3 &p0, const Vec3 &p1, const Vec3 &p2) {	
		v_min = min(min(p0, p1), p2);
        v_max = max(max(p0, p1), p2);
	}

    BoundingBox(const Vec3 &p0, const Vec3 &p1, const Vec3 &p2, const Vec3 &p3) {
		v_min = min(min(min(p0, p1), p2), p3);
        v_max = max(max(max(p0, p1), p2), p3);
	}
    
	/**
     * Return true if two boxes have intersection.
     */
	bool overlap(const BoundingBox *b) const {
        bool x = (v_max.x() >= b->v_min.x()) && (v_min.x() <= b->v_max.x());
        bool y = (v_max.y() >= b->v_min.y()) && (v_min.y() <= b->v_max.y());
        bool z = (v_max.z() >= b->v_min.z()) && (v_min.z() <= b->v_max.z());
        return (x && y && z);
        
        /*
        bool x = ((v_min.x() <= b->v_min.x() && b->v_min.x() <= v_max.x()) ||
                  (v_min.x() <= b->v_max.x() && b->v_max.x() <= v_max.x()) ||
                  (b->v_min.x() <= v_min.x() && v_max.x() <= b->v_max.x()));
        bool y = ((v_min.y() <= b->v_min.y() && b->v_min.y() <= v_max.y()) ||
                  (v_min.y() <= b->v_max.y() && b->v_max.y() <= v_max.y()) ||
                  (b->v_min.y() <= v_min.y() && v_max.y() <= b->v_max.y()));
        bool z = ((v_min.z() <= b->v_min.z() && b->v_min.z() <= v_max.z()) ||
                  (v_min.z() <= b->v_max.z() && b->v_max.z() <= v_max.z()) ||
                  (b->v_min.z() <= v_min.z() && v_max.z() <= b->v_max.z()));

        return (x && y && z);
        */
    }

	void merge(const BoundingBox &b) {
		v_min = min(v_min, b.v_min);
		v_max = max(v_max, b.v_max);		
	}

    void merge(Vec3 &p) {
		v_min = min(v_min, p);
		v_max = max(v_max, p);		
	}

    friend BoundingBox merge(const BoundingBox &a, const BoundingBox &b) {
        return BoundingBox(
            min(a.v_min, b.v_min),
		    max(a.v_max, b.v_max)
        );		
    }

    bool contains(const BoundingBox &b) const {
        bool x = ((v_min.x() <= b.v_min.x() && b.v_max.x() <= v_max.x()));
        bool y = ((v_min.y() <= b.v_min.y() && b.v_max.y() <= v_max.y()));
        bool z = ((v_min.z() <= b.v_min.z() && b.v_max.z() <= v_max.z()));

        return (x && y && z);
    }

    bool contains(const Vec3 &v) const {
        bool x = ((v_min.x() <= v.x() && v.x() <= v_max.x()));
        bool y = ((v_min.y() <= v.y() && v.y() <= v_max.y()));
        bool z = ((v_min.z() <= v.z() && v.z() <= v_max.z()));

        return (x && y && z);
    }

    /**
     * Scale v into range [-1, 1]
     */
    Vec3 normalize(const Vec3 &v) const {
        return (2.0f * v - v_max - v_min) / diagonal();
    }

    void reset() {
        v_min = Vec3(0.f);
        v_max = Vec3(0.f);
    }

    Vec3 size() const {
        return v_max - v_min;
    }

    Float diagonal() const {
        return (v_max - v_min).length();
    }

    Vec3 centroid() const {
        return 0.5f * (v_max + v_min);
    }

    bool hit(const Ray &r, Float tmin, Float tmax) const {
        Float t;
        return hit(r, tmin, tmax, t);
    }

    /*
    bool hit(const Ray &r, Float tmin, Float tmax, Float &t_hit) const {
        // transform so that the bounding box min is origin. 
        Vec3 p = r.org() - v_min;
        Vec3 d = r.dir();
        Vec3 size = v_max - v_min;

		// special case: when d is parallel to plane x = alpha, 
		// the ray can only hit the box if point p is inside the left and right slab.
		if (d.x() == 0) {
			if (p.x() < 0 || p.x() >= size.x()) return false;            
        }

		if (d.y() == 0)
			if (p.y() < 0 || p.y() >= size.y()) return false;

		if (d.z() == 0)
			if (p.z() < 0 || p.z() >= size.z()) return false;

        // intersection with pairs of plane x = 0 and x = size.x(), and so on.         
		// NOTE: might have division by zero here, but we will skip checking the corrupted t0 and t1
		// in the loop later.
        //Vec3 t0 = -p / d;
        //Vec3 t1 = (size - p) / d;  
        
        // Kay - Kajiya's method in Introduction to Ray Tracing by Andrew Glassner, pg. 65.
        Float t_x0, t_x1, t_y0, t_y1, t_z0, t_z1;
        Float t[6];
        int num = 0;

        if (d.x() != 0) {
            t_x0 = -p.x() / d.x();
            t_x1 = (size.x() - p.x()) / d.x();
            if (t_x0 > t_x1) std::swap(t_x0, t_x1);

            t[num++] = t_x0; // add to array to check
            t[num++] = t_x1;
        }

        if (d.y() != 0) {
            t_y0 = -p.y() / d.y();
            t_y1 = (size.y() - p.y()) / d.y();
            if (t_y0 > t_y1) std::swap(t_y0, t_y1);

            t[num++] = t_y0;
            t[num++] = t_y1;
        }

        if (d.z() != 0) {
            t_z0 = -p.z() / d.z();
            t_z1 = (size.z() - p.z()) / d.z();
            if (t_z0 > t_z1) std::swap(t_z0, t_z1);

            t[num++] = t_z0;
            t[num++] = t_z1;
        }

        if (num == 0) { // it seems ray direction is all zero 
            return false;
        }

        Float t_near = -FLT_MAX;
        Float t_far  =  FLT_MAX;        
        for (int i = 0; i < num; i += 2) {
            Float t1 = t[i    ];
            Float t2 = t[i + 1];

            if (t1 > t_near) t_near = t1;
            if (t2 < t_far) t_far = t2;

            if (t_near > t_far) return false; // box missed
            if (t_far < 0) return false; // box behind the ray
        }
        t_hit = t_near;
        if (t_hit < tmin) t_hit = t_far;
        
        //return (t_hit >= tmin && t_hit <= tmax); // t_hit <= tmax will think hit will be false when the segment is in the box.
        return (t_hit >= tmin); // as long as the ray segment hits the box, or inside the box, further investigation into the box is needed. Therefore, cannot bound tmax.
    }
    */
    
    /**
     * Find horizontal min and max in a __m128 vector.
     */
    /*
    float hor_min(__m128 val) const {
        // [A, B, C, D (index 0, 1, 2, 3) --> (C, D, A, B) -> min -> (min(A, C), min(B, D), min (A, C), min(B, D))
        // -> shuffle one more time and min one more time
        __m128 tmp = _mm_shuffle_ps(val, val, _MM_SHUFFLE(2, 3, 0, 1));
        val = _mm_min_ps (val, tmp);    
        tmp = _mm_shuffle_ps(val, val, _MM_SHUFFLE(1, 0, 3, 2));
        tmp = _mm_min_ps (val, tmp);
        return tmp.m128_f32[3];
    }
    float hor_max(__m128 val) const {        
        __m128 tmp = _mm_shuffle_ps(val, val, _MM_SHUFFLE(2, 3, 0, 1));
        val = _mm_max_ps (val, tmp);    
        tmp = _mm_shuffle_ps(val, val, _MM_SHUFFLE(1, 0, 3, 2));
        tmp = _mm_max_ps (val, tmp);
        return tmp.m128_f32[3];
    }
    */
    /*
    float hor_min(__m128 val) const {
        return std::min(val.m128_f32[0],
               std::min(val.m128_f32[1],
               std::min(val.m128_f32[2], 
                        val.m128_f32[3])));
    }
    float hor_max(__m128 val) const {        
        return std::max(val.m128_f32[0],
               std::max(val.m128_f32[1],
               std::max(val.m128_f32[2], 
                        val.m128_f32[3])));
    }*/

    bool hit(const Ray &r, Float tmin, Float tmax, Float &t_hit) const {
        // NOTE: how to handle NaN or Inf?
        // we handle it by setting inv_d to a large value when a ray is created. 

        // precompute for each ray
        //Vec3 inv_d = Vec3(1.0f, 1.0f, 1.0f) / r.dir();
        
        // from Timo Aila's HPG'12 Understanding the efficiency of ray traversal on the GPUs
        /*
        Vec3 v0 = (v_min - r.org()) * r.inv_d;
        Vec3 v1 = (v_max - r.org()) * r.inv_d;
        Float x0 = v0.x(), y0 = v0.y(), z0 = v0.z(),
              x1 = v1.x(), y1 = v1.y(), z1 = v1.z();
        */
        Vec3 v0 = v_min * r.inv_d + r.org_inv_d;  
        Vec3 v1 = v_max * r.inv_d + r.org_inv_d;
        //Vec3 v0 = mul_add(v_min, r.inv_d, r.org_inv_d);    // mul_add instruction not available
        //Vec3 v1 = mul_add(v_max, r.inv_d, r.org_inv_d);
        Float x0 = v0.x(), y0 = v0.y(), z0 = v0.z(),
              x1 = v1.x(), y1 = v1.y(), z1 = v1.z();
                
        //using (a < b ? a : b) instead of std::min is slower.                
        #define rz_min2(a, b) std::min(a, b)
        #define rz_max2(a, b) std::max(a, b)        
        #define rz_min4(a, b, c, d) (rz_min2(a, rz_min2(b, rz_min2(c, d))))
        #define rz_max4(a, b, c, d) (rz_max2(a, rz_max2(b, rz_max2(c, d))))
        Float tminbox = rz_max4(tmin, rz_min2(x0, x1), rz_min2(y0, y1), rz_min2(z0, z1));
        Float tmaxbox = rz_min4(tmax, rz_max2(x0, x1), rz_max2(y0, y1), rz_max2(z0, z1));
        
        /*
        // SSE and horizontal min/max is slower        
        __m128 m0 = _mm_load_ps(v0.e);
        __m128 m1 = _mm_load_ps(v1.e);
        __m128 smaller = _mm_min_ps(m0, m1); smaller.m128_f32[3] = tmin;
        __m128 larger = _mm_max_ps(m0, m1); larger.m128_f32[3] = tmax;
        Float tminbox = hor_max(smaller);
        Float tmaxbox = hor_min(larger);
        */

        t_hit = tminbox;
        return (tminbox <= tmaxbox) && 
               (tmax >= tminbox && tmaxbox >= tmin); /* some overlap of [tmin, tmax] and [tminbox, tmaxbox] */
    }
    	
    /**
     * Find minimum squared distance between point P and the box.
     */
    Float min_squared_distance(const Vec3 &p) const {
        Vec3 clamped_p = clamp(p, v_min, v_max);
        return (clamped_p - p).squared_length();
    }
        
    Float min_squared_distance1(const Vec3 &p) const {
        // TODO: from lightslice code, 
        // this is simply a clamp operation to apply on p with [v_min, v_max].
        // The minimum distance is the distance between p and the clamped p.

        // TODO: a faster approach is to check if p falls into the region as the diagram below
        // and determine the distance. The region are mutual exclusive.
        // B | 3 | C
        // ----------
        // 2 |   | 4
        // ----------
        // A | 1 | D

        Float min_dist = FLT_MAX;
        Float dist;
        if (contains(p)) {
            // test with six faces 
            dist = (p.x() - v_min.x()) * (p.x() - v_min.x());
            if (dist < min_dist)
                min_dist = dist;

            dist = (p.y() - v_min.y()) * (p.y() - v_min.y());
            if (dist < min_dist)
                min_dist = dist;

            dist = (p.z() - v_min.z()) * (p.z() - v_min.z());
            if (dist < min_dist)
                min_dist = dist;

            dist = (p.x() - v_max.x()) * (p.x() - v_max.x());
            if (dist < min_dist)
                min_dist = dist;

            dist = (p.y() - v_max.y()) * (p.y() - v_max.y());
            if (dist < min_dist)
                min_dist = dist;

            dist = (p.z() - v_max.z()) * (p.z() - v_max.z());
            if (dist < min_dist)
                min_dist = dist;

            return min_dist;
        }

        // if P does not fall into the faces of the box, check only the eight corners.
        // otherwise, check the face projections.

        // v_min planes
        if (contains(Vec3(v_min.x(), p.y(), p.z()))) {
            min_dist = p.x() - v_min.x();
            return min_dist * min_dist;
        } 
        if (contains(Vec3(p.x(), v_min.y(), p.z()))) {
            min_dist = p.y() - v_min.y();
            return min_dist * min_dist;
        }
        if (contains(Vec3(p.x(), p.y(), v_min.z()))) {
            min_dist = p.z() - v_min.z();
            return min_dist * min_dist;
        }
        // v_max planes
        if (contains(Vec3(v_max.x(), p.y(), p.z()))) {
            min_dist = p.x() - v_max.x();
            return min_dist * min_dist;
        } 
        if (contains(Vec3(p.x(), v_max.y(), p.z()))) {
            min_dist = p.y() - v_max.y();
            return min_dist * min_dist;
        }
        if (contains(Vec3(p.x(), p.y(), v_max.z()))) {
            min_dist = p.z() - v_max.z();
            return min_dist * min_dist;
        }
        
        // check corners
        dist = (p - v_min).squared_length();
        if (dist < min_dist)
            min_dist = dist;

        dist = (p - Vec3(v_min.x(), v_min.y(), v_max.z())).squared_length();
        if (dist < min_dist)
            min_dist = dist;

        dist = (p - Vec3(v_min.x(), v_max.y(), v_min.z())).squared_length();
        if (dist < min_dist)
            min_dist = dist;
     
        dist = (p - Vec3(v_min.x(), v_max.y(), v_max.z())).squared_length();
        if (dist < min_dist)
            min_dist = dist;

        dist = (p - Vec3(v_max.x(), v_min.y(), v_min.z())).squared_length();
        if (dist < min_dist)
            min_dist = dist;

        dist = (p - Vec3(v_max.x(), v_min.y(), v_max.z())).squared_length();
        if (dist < min_dist)
            min_dist = dist;

        dist = (p - Vec3(v_max.x(), v_max.y(), v_min.z())).squared_length();
        if (dist < min_dist)
            min_dist = dist;

        dist = (p - v_max).squared_length();
        if (dist < min_dist)
            min_dist = dist;            

        return min_dist;  
    }

public:
	Vec3 v_min, v_max;
};

typedef vector<BoundingBox> BoundingBoxes;


class BoundingSphere {
public:
    BoundingSphere() : rad(-1.0f) {}
    
    void reset() {
        rad = -1.0f;
    }

    BoundingSphere(const Vec3 &center, Float rad) : center(center), rad(rad) {}
    
    BoundingSphere(const Vec3 &x, const Vec3 &y) {
        rad = 0.5f * (x - y).length();
        center = 0.5f * (x + y);
    }

    BoundingSphere(const Vec3 &a, const Vec3 &b, const Vec3 &c) {
        Onb uvn;
        uvn.init_from_uv(unit_vector(b - a), unit_vector(c - a));
        Vec3 aa = uvn.world_to_local(a);
        Vec3 bb = uvn.world_to_local(b);
        Vec3 cc = uvn.world_to_local(c);
        
        // translate to A
        bb -= aa;
        cc -= aa;

        // find circumcenter
        Float dd = 2.0f * (bb.x() * cc.y() - bb.y() * cc.x());
        if (dd == 0.0f) {
            // degenerate triangle, three points colinear
            rad = 0.5f * (a - b).length();
            center = 0.5f * (a + b);
            this->merge(c);
            return;
        }

        // equation from Wikipedia: http://en.wikipedia.org/wiki/Circumscribed_circle
        Float center_x = (cc.y() * (bb.x() * bb.x() + bb.y() * bb.y()) - bb.y() * (cc.x() * cc.x() + cc.y() * cc.y())) / dd;
        Float center_y = (bb.x() * (cc.x() * cc.x() + cc.y() * cc.y()) - cc.x() * (bb.x() * bb.x() + bb.y() * bb.y())) / dd;
        Vec3 circumcenter(center_x, center_y, 0.0f);
        
        center = uvn.local_to_world(circumcenter);        
        rad = sqrt(center_x * center_x + center_y * center_y);          // because A is at origin.

        // test
        Float ra = (center - a).length();
        Float rb = (center - b).length();
        Float rc = (center - c).length();
        if (fabs(ra - rb) > ZERO_EPSILON || 
            fabs(ra - rc) > ZERO_EPSILON ||
            fabs(rb - rc) > ZERO_EPSILON) {
                Log::info() << "Error: "  << ra << " " << rb << " " << rc << endn;
        }
    }

    //bool hit(const Ray &r, Float tmin, Float tmax, Float time) const;

    void merge(const Vec3 &p) {
        if (rad < 0.0f) {
            center = p;
            rad = 0.0f;
            return;
        }

        if (contains(p)) return;

        rad = 0.5f * (p - center).length();
        center = p + unit_vector(center - p) * rad;
    }

    void merge(const BoundingSphere &s) {
        if (rad < 0.0f) {
            center = s.center;
            rad = s.rad;
            return;
        }

        if (contains(s)) return;

        if (s.contains(*this)) {
            rad = s.rad;
            center = s.center;
            return;
        }

        rad = 0.5f * ((s.center - center).length() + s.rad + rad);
        center = s.center + (rad - s.rad) * (center - s.center);
    }

    bool overlap(const BoundingSphere &s) const {
        if (rad < 0.0f) return false;

        return (center - s.center).length() < rad + s.rad;
    }

    bool contains(const Vec3 &p) const {
        if (rad < 0.0f) return false;

        return ((p - center).length() <= rad);
    }

    bool contains(const BoundingSphere &s) const {
        if (rad < 0.0f) return false;

        Vec3 test_point = s.center + unit_vector(s.center - center) * s.rad;
        return contains(test_point);
    }

public:
    Vec3 center;
    Float rad;
};


} // end namespace

#endif
