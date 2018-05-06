#include "triangle.h"

namespace Renzoku {

Triangle::Triangle() {}

void Triangle::init(const Vec3 &_p0, const Vec3 &_p1, const Vec3 &_p2) {    
    p0 = _p0;
    p1 = _p1;
    p2 = _p2;

    Vec3 p01 = p0 - p1;
    Vec3 p02 = p0 - p2;
    
    tri_area = cross(p01, p02).length() * 0.5f;    
    
    u = unit_vector(p01);
    n = unit_vector(cross(p01, p02));
    n0 = n1 = n2 = n;

    compute_bounding_box();
}

Triangle::Triangle(const Vec3 &_p0, const Vec3 &_p1, const Vec3 &_p2) {    
    init(_p0, _p1, _p2);
}

Triangle::Triangle(const Vec3 &_p0, const Vec3 &_p1, const Vec3 &_p2,
                   const Vec2 &uv0, const Vec2 &uv1, const Vec2 &uv2) 
      : uv0(uv0), uv1(uv1), uv2(uv2)
{
    init(_p0, _p1, _p2);
}

/**
 * Follow algorithm by Didier Badouel at:
 * http://graphics.stanford.edu/courses/cs348b-98/gg/intersect.html
 */
 /*
bool Triangle::hit(const Ray &r, Float tmin, Float tmax, Float time, HitRecord &record) const {
    Vec3 v01 = p1 - p0;
    Vec3 v02 = p2 - p0;

    // plane of triangle: n . p + d = 0
    Vec3  n = cross(v01, v02);
    Float d = - dot(n, p0);

    Vec3 O = r.org();
    Vec3 D = r.dir();
    Float dot_n_D = dot(n, D);
    if (dot_n_D == 0) return false;
    
    Float t = - (d + dot(n, O)) / dot_n_D;
    if (t < tmin || t > tmax) return false;
    
    // find dominant axis of the normal
    Float abs_nx = fabs(n.x());
    Float abs_ny = fabs(n.y());
    Float abs_nz = fabs(n.z());
    int i1, i2;
    
    // WRONG: 
    if (abs_nx > abs_ny) {
        i1 = 1;
        if (abs_nx > abs_nz) { // x max            
            i2 = 2;
        } else { // z max
            i2 = 0;            
        }
    } else {
        i1 = 0;
        if (abs_ny > abs_nz) { // y max            
            i2 = 2;
        } else { // z max
            i2 = 1;
        }
    }

    Vec3 p = O + t * D;
    Float u0 = p[i1] - p0[i1];
    Float v0 = p[i2] - p0[i2];

    Float u1 = p1[i1] - p0[i1];
    Float v1 = p1[i2] - p0[i2];

    Float u2 = p2[i1] - p0[i1];
    Float v2 = p2[i2] - p0[i2];

    Float det = u1 * v2 - u2 * v1;
    Float alpha = (u0 * v2 - u2 * v0) / det;
    Float beta  = (u1 * v0 - u0 * v1) / det;
    
    bool hit = (alpha >= 0 && beta >= 0 && (alpha + beta <= 1)); 
    if (hit) {
        record.t = t;
        record.normal = unit_vector(n);
        //cout << "Tri: " << p01 << " " << p02 << " " << record.normal << endl;
        record.uvn->init_from_nu(n, v01);
    }
    return hit;
}

bool Triangle::hit(const Ray &r, Float tmin, Float tmax, Float time) const {
    Vec3 v01 = p1 - p0;
    Vec3 v02 = p2 - p0;

    // plane of triangle: n . p + d = 0
    Vec3  n = cross(v01, v02);
    Float d = - dot(n, p0);

    Vec3 O = r.org();
    Vec3 D = r.dir();
    Float dot_n_D = dot(n, D);
    if (dot_n_D == 0) return false;
    
    Float t = - (d + dot(n, O)) / dot_n_D;
    if (t < tmin || t > tmax) return false;
    
    // find dominant axis of the normal
    Float abs_nx = fabs(n.x());
    Float abs_ny = fabs(n.y());
    Float abs_nz = fabs(n.z());
    int i1, i2;
    
    if (abs_nx > abs_ny) {
        i1 = 1;
        if (abs_nx > abs_nz) { // x max            
            i2 = 2;
        } else { // z max
            i2 = 0;            
        }
    } else {
        i1 = 0;
        if (abs_ny > abs_nz) { // y max            
            i2 = 2;
        } else { // z max
            i2 = 1;
        }
    }

    Vec3 p = O + t * D;
    Float u0 = p[i1] - p0[i1];
    Float v0 = p[i2] - p0[i2];

    Float u1 = p1[i1] - p0[i1];
    Float v1 = p1[i2] - p0[i2];

    Float u2 = p2[i1] - p0[i1];
    Float v2 = p2[i2] - p0[i2];

    Float det = u1 * v2 - u2 * v1;
    Float alpha = (u0 * v2 - u2 * v0) / det;
    Float beta  = (u1 * v0 - u0 * v1) / det;
    
    return (alpha >= 0 && beta >= 0 && (alpha + beta <= 1));     
}*/

bool Triangle::hit(const Ray &r, Float tmin, Float tmax, Float time, GeometryHit &record) const {
    record.hit = false;

	/*
	 * Solve 3 linear equation p0 + t * d = alpha * AB + beta * AC
	 * to obtain t, alpha, beta.
     * 
     * Special case:
     * When the ray is parallel to the triangle, 
     * d is a linear combination of AB and AC, which causes the matrix [AB, AC, d] 
     * to be ill-conditioned and cannot be inverted. The determinant of this matrix is zero.
     * Solution becomes none if the ray is not coplanar with the triangle, or a line segment
     * if the ray is coplanar. 
     */
    Vec3 p01 = p0 - p1;
    Float A = p01.x();
    Float B = p01.y();
    Float C = p01.z();

    Vec3 p02 = p0 - p2;
    Float D = p02.x();
    Float E = p02.y();
    Float F = p02.z();

    Float G = r.dir().x();  
    Float H = r.dir().y();
    Float I = r.dir().z();
    
    Vec3 p0org = p0 - r.org();
    Float J = p0org.x();
    Float K = p0org.y();    
    Float L = p0org.z();
    
    Float EIHF = E * I - H * F;
    Float GFDI = G * F - D * I;
    Float DHEG = D * H - E * G;
    
    Float denom = (A * EIHF + B * GFDI + C * DHEG);
	if (-ZERO_EPSILON < denom && denom < ZERO_EPSILON) { // ray lies on the same surface with the triangle, matrix degenerated
		return false;

		/*
		// test if ray is coplanar
		Vec3 n = unit_vector(cross(p01, p02));
		if (fabs(dot(n, p0org)) < epsilon) { // coplanar
			record.t = t_near;
			record.normal = n;
			record.uvn->init_from_nu(record.normal, p01);
			return true;
		} else {
			return false;
		}*/
	}

    Float beta = (J * EIHF + K * GFDI + L * DHEG) / denom;
    if (beta < 0.0f || beta > 1.0f) return false;

    Float AKJB = A * K - J * B;
    Float JCAL = J * C - A * L;
    Float BLKC = B * L - K * C;
    
    Float gamma = (I * AKJB + H * JCAL + G * BLKC) / denom;
    if (gamma < 0.0f || beta + gamma > 1.0f) return false;

    Float t = -(F * AKJB + E * JCAL + D * BLKC) / denom;
    if (t >= tmin && t <= tmax) {
        record.hit = true;
        record.t = t;        
        record.beta = beta;
        record.gamma = gamma;
        record.shape = (Shape *)this;
        return true;
    }   
    return false;
}

bool Triangle::hit(const Ray &r, Float tmin, Float tmax, Float time) const {
    Vec3 p01 = p0 - p1;
    Float A = p01.x();
    Float B = p01.y();
    Float C = p01.z();

    Vec3 p02 = p0 - p2;
    Float D = p02.x();
    Float E = p02.y();
    Float F = p02.z();

    Float G = r.dir().x();  
    Float H = r.dir().y();
    Float I = r.dir().z();
    
    Vec3 p0org = p0 - r.org();
    Float J = p0org.x();
    Float K = p0org.y();    
    Float L = p0org.z();
    
    Float EIHF = E * I - H * F;
    Float GFDI = G * F - D * I;
    Float DHEG = D * H - E * G;
    
    Float denom = (A * EIHF + B * GFDI + C * DHEG);
	if (-ZERO_EPSILON < denom && denom < ZERO_EPSILON) {
		return false;
	}

    Float beta = (J * EIHF + K * GFDI + L * DHEG) / denom;
    if (beta < 0.0f || beta > 1.0f) return false;

    Float AKJB = A * K - J * B;
    Float JCAL = J * C - A * L;
    Float BLKC = B * L - K * C;
    
    Float gamma = (I * AKJB + H * JCAL + G * BLKC) / denom;
    if (gamma < 0.0f || beta + gamma > 1.0f) return false;

    Float t = -(F * AKJB + E * JCAL + D * BLKC) / denom;
    return (t >= tmin && t <= tmax);
}

void Triangle::fill_hit_record(const Ray &r, const GeometryHit &gh, HitRecord &record) const {
    record.copy_geometry_hit(gh);

    record.normal         = this->n;
        
    // default shading normal is interpolated from vertex normals
    record.shading_normal = unit_vector(n0 + record.beta * (n1 - n0) + record.gamma * (n2 - n0));
                
    /**
     * TODO: 
     * 1. tangent must derive from texture UV for correct orientation. 
     */
    record.tangent = this->u;
        
    // beta for (p1 - p0)
    // gamma for (p2 - p0)
    record.uv = uv0 + record.beta * (uv1 - uv0) + record.gamma * (uv2 - uv0);

    //record.shape = (Shape *)this;
}

void Triangle::sample(Random &rd, Vec3 &p, Vec3 &n, Float &pdf) const {
    Vec2 uv;
    uv.random(rd);
    Float tmp = sqrt(1.0f - uv[0]);
    Float b = 1.0f - tmp;
    Float c = uv[1] * tmp;
    Float a = 1.0f - b - c;
    p = (a * p0 + b * p1 + c * p2);
    n = this->normal(p);
    pdf = 1.0f / tri_area;
}

Float Triangle::pdf(const Vec3 &p) const {
    Float a, b, c;
    Linear::solve33(p0.x(), p1.x(), p2.x(),
                    p0.y(), p1.y(), p2.y(),
                    p0.z(), p1.z(), p2.z(),
                    p.x(), p.y(), p.z(), 
                    a, b, c);

    if (a >= 0.0f && a <= 1.0f && 
        b >= 0.0f && b <= 1.0f &&
        c >= 0.0f && c <= 1.0f) {

        return 1.0f / tri_area;

    } else {
        return 0.0f;
    }
}

void Triangle::get_triangle_coordinates(vector<Float> &coords) const {
    coords.push_back(p0.x());
    coords.push_back(p0.y());
    coords.push_back(p0.z());
    coords.push_back(p1.x());
    coords.push_back(p1.y());
    coords.push_back(p1.z());
    coords.push_back(p2.x());
    coords.push_back(p2.y());
    coords.push_back(p2.z());
}

void Triangle::get_vertex_data(vector<GlVertex> &vertices) const {    
    GlVertex v;    
    v.normal[0] = n0.x();
    v.normal[1] = n0.y();
    v.normal[2] = n0.z();
    v.pos[0] = p0.x();
    v.pos[1] = p0.y();
    v.pos[2] = p0.z();
    v.tangent[0] = u.x();
    v.tangent[1] = u.y();
    v.tangent[2] = u.z();
    v.uv[0] = uv0.x();
    v.uv[1] = uv0.y();
    vertices.push_back(v);

    v.normal[0] = n1.x();
    v.normal[1] = n1.y();
    v.normal[2] = n1.z();
    v.pos[0] = p1.x();
    v.pos[1] = p1.y();
    v.pos[2] = p1.z();
    v.tangent[0] = u.x();
    v.tangent[1] = u.y();
    v.tangent[2] = u.z();
    v.uv[0] = uv1.x();
    v.uv[1] = uv1.y();
    vertices.push_back(v);

    v.normal[0] = n2.x();
    v.normal[1] = n2.y();
    v.normal[2] = n2.z();
    v.pos[0] = p2.x();
    v.pos[1] = p2.y();
    v.pos[2] = p2.z();
    v.tangent[0] = u.x();
    v.tangent[1] = u.y();
    v.tangent[2] = u.z();
    v.uv[0] = uv2.x();
    v.uv[1] = uv2.y();
    vertices.push_back(v);
}

void Triangle::get_vertex_positions(vector<Vec3> &positions) const {
    positions.push_back(p0);
    positions.push_back(p1);
    positions.push_back(p2);
}

void Triangle::set_vertex_normal(int vertex_idx, const Vec3 &normal) {
    switch (vertex_idx) {
    case 0: 
        n0 = normal;
        break;
    case 1:
        n1 = normal;
        break;
    case 2:
        n2 = normal;
        break;
    default:
        break;
    }
}

void Triangle::set_tangent(const Vec3 &tangent) {
    u = unit_vector(tangent);
}

} // end namespace Renzoku

