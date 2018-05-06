#include "cone.h"
#include "log.h"

namespace Renzoku {

Float Cone::angle(const Vec3 &a, const Vec3 &b) {
    // go from a to b anti-clockwise
    // assume 2D vector (x, y, 0)
  
    Float alpha = fmod(atan2(-b.x() * a.y() + a.x() * b.y(), a.x() * b.x() + a.y() * b.y()) + TWO_PI, TWO_PI);
    if (alpha != alpha) {
        Log::info() << "NaN angle 1" << endn;
    }
    return alpha;
}

Float Cone::angle_from_x_axis(const Vec3 &a) {
    // go from (1, 0) to a anti-clockwise
    // assume 2D vector (x, y, 0)
    
    Float alpha = fmod(atan2(a.y(), a.x()) + TWO_PI, TWO_PI);
    if (alpha != alpha) {
        Log::info() << "NaN angle 2" << endn;
    }
    return alpha;
}

Cone Cone::merge_direction(const Cone &src, const Cone &dest) {
    if ((src.normal() + dest.normal()).squared_length() == 0.0f &&
        src.half_angle() == 0.0f && dest.half_angle() == 0.0f) {

        // degenerate case, the bounding cone can rotate around the src and dest normal "axis"
        // we need to find the half vector by rotating the src normal
        Onb uvn;
        uvn.init_from_u(src.normal());    
        Vec3 n = uvn.local_to_world(rotate_about_z(uvn.world_to_local(src.normal()), HALF_PI));
        return Cone(n, HALF_PI);
    }

	Onb uvn;
	uvn.init_from_uv(src.n, dest.n);	
	Vec3 a_n = uvn.world_to_local(src.n);
    Vec3 b_n = uvn.world_to_local(dest.n);    
    // the vectors are now in the form (x, y, 0). 
    // sin of (x1, y1) to (x2, y2) is x2 * y1 - x1 * y2.

    //cout << "a_n: " << a_n << endl;

	Vec3 a1 = rotate_about_z(a_n, -src.half);
	Vec3 a2 = rotate_about_z(a_n,  src.half);
	Vec3 b1 = rotate_about_z(b_n, -dest.half);
	Vec3 b2 = rotate_about_z(b_n,  dest.half);
    
    Float aa1 = angle_from_x_axis(a1);
    Float aa2 = angle_from_x_axis(a2);
    Float bb1 = angle_from_x_axis(b1);
    Float bb2 = angle_from_x_axis(b2);

    // Assumption: aa1 < aa2, bb1 < bb2
    if (aa1 > aa2) {
        aa1 -= TWO_PI;
    }
    if (bb1 > bb2) {
        bb1 -= TWO_PI;
    }

    //cout << "aa1: " << aa1 * 180 / A_PI << endl;
    //cout << "aa2: " << aa2 * 180 / A_PI << endl;
    //cout << "bb1: " << bb1 * 180 / A_PI << endl;
    //cout << "bb2: " << bb2 * 180 / A_PI << endl;

    Vec3 n;
    Float half;

    if (bb1 <= aa1) {
        if (bb2 < aa1) {
            Float angle1 = angle(b1, a2);
            Float angle2 = angle(a1, b2);
            if (angle1 < angle2) {
                half = 0.5f * angle1;
                n = uvn.local_to_world(unit_vector((b1 + a2) * 0.5f));
                if (half == 0 && ((aa1 != aa2) || (bb1 != bb2))) half = A_PI;
                if (half > HALF_PI) n = -n;
            } else {
                half = 0.5f * angle2;
                n = uvn.local_to_world(unit_vector((a1 + b2) * 0.5f));
                if (half == 0 && ((aa1 != aa2) || (bb1 != bb2))) half = A_PI;
                if (half > HALF_PI) n = -n;
            }
        } else if (bb2 <= aa2) {
            half = angle(b1, a2) * 0.5f;
            n = uvn.local_to_world(unit_vector((b1 + a2) * 0.5f));
            if (half == 0 && ((aa1 != aa2) || (bb1 != bb2))) half = A_PI;
            if (half > HALF_PI) n = -n;
        } else {
            half = dest.half;
            n = dest.n;
        }
    } else if (bb1 <= aa2) {
        if (bb2 <= aa2) {
            half = src.half;
            n = src.n;
        } else {
            half = 0.5f * angle(a1, b2);
            n = uvn.local_to_world(unit_vector((a1 + b2) * 0.5f));
            if (half == 0 && ((aa1 != aa2) || (bb1 != bb2))) half = A_PI;
            if (half > HALF_PI) n = -n;
        }
    } else {
        Float angle1 = angle(b1, a2);
        Float angle2 = angle(a1, b2);
        //cout << "angle1: " << angle1 * 180 / A_PI << endl;
        //cout << "angle2: " << angle2 * 180 / A_PI << endl;
        if (angle1 < angle2) {
            half = 0.5f * angle1;
            n = uvn.local_to_world(unit_vector((b1 + a2) * 0.5f));
            if (half == 0 && ((aa1 != aa2) || (bb1 != bb2))) half = A_PI;
            if (half > HALF_PI) n = -n;
        } else {
            half = 0.5f * angle2;
            n = uvn.local_to_world(unit_vector((a1 + b2) * 0.5f));
            if (half == 0 && ((aa1 != aa2) || (bb1 != bb2))) half = A_PI;
            if (half > HALF_PI) n = -n;
        }
    }
        
    return Cone(n, half);
}

void Cone::merge_direction(const Cone &src) {
    if (this->n.squared_length() == 0.0f) {
        // when the cone is not initialized
        this->n = src.n;
        this->half = src.half;
        return;
    }

    if ((src.normal() + this->normal()).squared_length() == 0.0f &&
        src.half_angle() == 0.0f && this->half_angle() == 0.0f) {

        // degenerate case, the bounding cone can rotate around the src and dest normal "axis"
        // we need to find the half vector by rotating the src normal
        Onb uvn;
        uvn.init_from_u(src.normal());
        Vec3 n = uvn.local_to_world(rotate_about_z(uvn.world_to_local(src.normal()), HALF_PI));
       
        this->n = n;
        this->half = HALF_PI;
        return;
    }

    Onb uvn;
	uvn.init_from_uv(src.n, this->n);	
	Vec3 a_n = uvn.world_to_local(src.n);
    Vec3 b_n = uvn.world_to_local(this->n);    
    // the vectors are now in the form (x, y, 0). 
    // sin of (x1, y1) to (x2, y2) is x2 * y1 - x1 * y2.
    
	Vec3 a1 = rotate_about_z(a_n, -src.half);
	Vec3 a2 = rotate_about_z(a_n,  src.half);
	Vec3 b1 = rotate_about_z(b_n, -this->half);
	Vec3 b2 = rotate_about_z(b_n,  this->half);
    
    Float aa1 = angle_from_x_axis(a1);
    Float aa2 = angle_from_x_axis(a2);
    Float bb1 = angle_from_x_axis(b1);
    Float bb2 = angle_from_x_axis(b2);

    // Assumption: aa1 < aa2, bb1 < bb2
    if (aa1 > aa2) {
        aa1 -= TWO_PI;
    }
    if (bb1 > bb2) {
        bb1 -= TWO_PI;
    }
    
    Vec3 n;
    Float half;

    if (bb1 <= aa1) {
        if (bb2 < aa1) {
            Float angle1 = angle(b1, a2);
            Float angle2 = angle(a1, b2);
            if (angle1 < angle2) {
                half = 0.5f * angle1;
                n = uvn.local_to_world(unit_vector((b1 + a2) * 0.5f));
                if (half == 0 && ((aa1 != aa2) || (bb1 != bb2))) half = A_PI;
                if (half > HALF_PI) n = -n;
            } else {
                half = 0.5f * angle2;
                n = uvn.local_to_world(unit_vector((a1 + b2) * 0.5f));
                if (half == 0 && ((aa1 != aa2) || (bb1 != bb2))) half = A_PI;
                if (half > HALF_PI) n = -n;
            }
        } else if (bb2 <= aa2) {
            half = angle(b1, a2) * 0.5f;
            n = uvn.local_to_world(unit_vector((b1 + a2) * 0.5f));
            if (half == 0 && ((aa1 != aa2) || (bb1 != bb2))) half = A_PI;
            if (half > HALF_PI) n = -n;
        } else {
            // keep as is
            n = this->n;
            half = this->half;
        }
    } else if (bb1 <= aa2) {
        if (bb2 <= aa2) {
            half = src.half;
            n = src.n;
        } else {
            half = 0.5f * angle(a1, b2);
            n = uvn.local_to_world(unit_vector((a1 + b2) * 0.5f));
            if (half == 0 && ((aa1 != aa2) || (bb1 != bb2))) half = A_PI;
            if (half > HALF_PI) n = -n;
        }
    } else {
        Float angle1 = angle(b1, a2);
        Float angle2 = angle(a1, b2);
        
        if (angle1 < angle2) {
            half = 0.5f * angle1;
            n = uvn.local_to_world(unit_vector((b1 + a2) * 0.5f));
            if (half == 0 && ((aa1 != aa2) || (bb1 != bb2))) half = A_PI;
            if (half > HALF_PI) n = -n;
        } else {
            half = 0.5f * angle2;
            n = uvn.local_to_world(unit_vector((a1 + b2) * 0.5f));
            if (half == 0 && ((aa1 != aa2) || (bb1 != bb2))) half = A_PI;
            if (half > HALF_PI) n = -n;
        }
    }

    this->n = n;
    this->half = half;
}

void Cone::merge_direction(const Vec3 &v) {
    if (this->n.squared_length() == 0.0f) {
        // when the cone is not initialized
        this->n = v;
        this->half = 0.0f;
        return;
    }

    if ((v + this->normal()).squared_length() == 0.0f &&
        this->half_angle() == 0.0f) {

        // degenerate case, the bounding cone can rotate around the src and dest normal "axis"
        // we need to find the half vector by rotating the src normal
        Onb uvn;
        uvn.init_from_u(v);
        Vec3 n = uvn.local_to_world(rotate_about_z(uvn.world_to_local(v), HALF_PI));
       
        // should we just keep the old n?

        this->n = n;
        this->half = HALF_PI;
        return;
    }
    
    if (contains(v)) return;

    Onb uvn;
	uvn.init_from_uv(v, this->n);		
    Vec3 a = uvn.world_to_local(v);
    Vec3 b_n = uvn.world_to_local(this->n);    
    // the vectors are now in the form (x, y, 0). 
    // sin of (x1, y1) to (x2, y2) is x2 * y1 - x1 * y2.    
	Vec3 b1 = rotate_about_z(b_n, -this->half);
	Vec3 b2 = rotate_about_z(b_n,  this->half);
    
    Float bb1 = angle_from_x_axis(b1);
    Float bb2 = angle_from_x_axis(b2);

    // Assumption: aa1 < aa2, bb1 < bb2
    if (bb1 > bb2) {
        bb1 -= TWO_PI;
    }

    Float angle1 = angle(b1, a);
    Float angle2 = angle(a, b2);
        
    if (angle1 < angle2) {
        half = 0.5f * angle1;
        n = uvn.local_to_world(unit_vector((b1 + a) * 0.5f));
        if (half == 0 && (bb1 != bb2)) half = A_PI;
        if (half > HALF_PI) n = -n;
    } else {
        half = 0.5f * angle2;
        n = uvn.local_to_world(unit_vector((a + b2) * 0.5f));
        if (half == 0 && ((bb1 != bb2))) half = A_PI;
        if (half > HALF_PI) n = -n;
    }
    this->half = half;
    this->n = n;
}

void Cone::merge_direction_algebraic(const Vec3 &v) {
    if (this->n.squared_length() == 0.0f) {
        // when the cone is not initialized
        this->n = v;
        this->half = 0.0f;
        return;
    }

    if ((v + this->normal()).squared_length() == 0.0f &&
        this->half_angle() == 0.0f) {

        // degenerate case, the bounding cone can rotate around the src and dest normal "axis"
        // we need to find the half vector by rotating the src normal
        Onb uvn;
        uvn.init_from_u(v);
        Vec3 n = uvn.local_to_world(rotate_about_z(uvn.world_to_local(v), HALF_PI));
       
        // should we just keep the old n?

        this->n = n;
        this->half = HALF_PI;
        return;
    }
    
    if (contains(v)) return;

    // we try to follow the approach by 
    // T.W. Sederberg and R.J. Meyers, Loop detection in surface patch intersections,
    // Computer-Aided Geometric Design, 5 (1998), 161-171.
    // 
    // the basic idea is to compute a bounding cone that contains the old cone and the new direction.     
    Vec3 nn = n;
    Float theta = half;
    /*
    if (half > HALF_PI) {
        nn = -n;
        theta = A_PI - half;
    }

    Cone opposite(-nn, theta);
    if (opposite.contains(v)) {
        nn = -nn;
        theta = A_PI - theta;
    }

    if (dot(v, nn) < 0.0f) {
        nn = -nn;
        theta = A_PI - theta;
    }*/

    Float cos_alpha = dot(v, nn);
    Float tan_theta = tan(theta);
    Vec3 Dt = unit_vector(nn + tan_theta * unit_vector(nn - v / cos_alpha));
    Vec3 new_n = unit_vector(Dt + v);
    Float new_half = acos(dot(new_n, v));

    if (new_half < half) {
        new_n = -new_n;
        new_half = acos(dot(new_n, v));
    }

    this->n = new_n;
    this->half = new_half;

    // this approach gives similar results to my version which actually implements this idea numerically.
    // but the equation now only works with nonreflex cones. So not using it now. 
}


} // end namespace