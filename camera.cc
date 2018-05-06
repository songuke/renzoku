#include "camera.h"
#include "scene.h"
#include "image.h"
#include "integrator.h"

namespace Renzoku {
Camera::Camera() : is_log_pixel(false), 
    near_plane(0.0f), far_plane(0.0f), focal(0.0f), scene(NULL),
    sensor(NULL) {
}

Camera::Camera(const Vec3 &_eye, const Vec3 &_lookat, const Vec3 &_up, Float _focal, Sensor *sensor)
    : eye(_eye), lookat(_lookat), up(_up), 
      focal(_focal),
      is_log_pixel(false),
      near_plane(0.0f), far_plane(0.0f)
{   
    onb.init_from_nv(eye - lookat, up); // looking at negative Z-axis

    this->set_sensor(sensor);
}

void Camera::set_perspective(const Vec3 &_eye, const Vec3 &_lookat, const Vec3 &_up, Float _focal, Sensor *sensor)
{    
    eye = _eye;
    lookat = _lookat;
    up = _up;
    focal = _focal;            
    onb.init_from_nv(eye - lookat, up); // looking at negative Z-axis
    this->set_sensor(sensor);
}

void Camera::set_perspective(const Vec3 &_eye, const Vec3 &_lookat, const Vec3 &_up, Float fov) {
    eye = _eye;
    lookat = _lookat;
    up = _up;
    onb.init_from_nv(eye - lookat, up); // looking at negative Z-axis

    // convert vertical field of view to focal
    Size2 film_size = sensor->get_film_size();
    focal = film_size.height / (2 * tan(fov / 2 * DEGREE_TO_RADIAN));    
    // NOTE: the aspect ratio together with vertical fov can infer the film size. 
}

Camera *Camera::clone() {
    Camera * c = new Camera(eye, lookat, up, focal, sensor);
    c->set_near_plane(near_plane);
    c->set_far_plane(far_plane);
    return c;
}

Ray Camera::shoot_ray(int i, int j) {
    Size2 img_size = sensor->get_image_size();
    Vec2 cop = sensor->get_cop();
    Vec2 scale = sensor->get_image_film_ratio();
    Vec2 pixel(j + 0.5f, img_size.height - 1 - i + 0.5f);
    Vec2 film_pixel = (pixel - cop) / scale;
    
    // direction in camera space since eye is origin
    Vec3 dir(film_pixel.x(), film_pixel.y(), -focal);  

    // direction in world space
    dir = unit_vector(onb.local_to_world(dir));    
    return Ray(eye, dir);
}

Ray Camera::shoot_ray(const SensorSample &sample) {
    Size2 img_size = sensor->get_image_size();
    Vec2 cop = sensor->get_cop();
    Vec2 scale = sensor->get_image_film_ratio();
    Vec2 pixel(sample.x + sample.jitter.x(), img_size.height - (sample.y + sample.jitter.y()));
    Vec2 film_pixel = (pixel - cop) / scale;
    
    // direction in camera space since eye is origin
    Vec3 dir(film_pixel.x(), film_pixel.y(), -focal);  

    // direction in world space
    dir = unit_vector(onb.local_to_world(dir));    
    return Ray(eye, dir);
}

Ray Camera::shoot_ray(int i, int j, Random &rd) {    
    Size2 img_size = sensor->get_image_size();
    Vec2 cop = sensor->get_cop();
    Vec2 scale = sensor->get_image_film_ratio();
    // jitter sampling
    Vec2 jitter; 
    jitter.random(rd);
    Vec2 pixel(j + jitter.x(), img_size.height - 1 - i + jitter.y());
    Vec2 film_pixel = (pixel - cop) / scale;
    
    // direction in camera space since eye is origin
    Vec3 dir(film_pixel.x(), film_pixel.y(), -focal);  

    // direction in world space
    dir = unit_vector(onb.local_to_world(dir));    
    return Ray(eye, dir);
}

Ray Camera::shoot_ray(int i, int j, Vec2 jitter) {
    Size2 img_size = sensor->get_image_size();
    Vec2 cop = sensor->get_cop();
    Vec2 scale = sensor->get_image_film_ratio();
    Vec2 pixel(j + jitter.x(), img_size.height - 1 - i + jitter.y());
    Vec2 film_pixel = (pixel - cop) / scale;
    
    // direction in camera space since eye is origin
    Vec3 dir(film_pixel.x(), film_pixel.y(), -focal);  

    // direction in world space
    dir = unit_vector(onb.local_to_world(dir));    
    return Ray(eye, dir);
}

bool Camera::hit(const Vec3 &p, Vec2 &pixel, Vec3 &wo_e) 
{
	// find the intersection point between the line(eye, p) and the film plane.    
    Vec3 u = onb.u();
    Vec3 v = onb.v();
    Vec3 n = onb.n();    

    // point p in camera coordinates
    Vec3 p0(dot(u, p - eye), dot(v, p - eye), dot(n, p - eye));

    // projected point: P' = E + t (P - E) where P is the 3D point; E is eye.
    // in camera coordinates, E = 0 --> P' = t * P    
    Float t = (-focal) / p0.z();
    Float xp = t * p0.x();
    Float yp = t * p0.y();
    
    // find the corresponding image pixel
    Vec2 filmPixel2(xp, yp);
    Vec2 cop = sensor->get_cop();
    Vec2 scale = sensor->get_image_film_ratio();
    pixel = filmPixel2 * scale + cop;
    
    // NOTE:    
    // The following code does the same thing, but lose precision in 32-bit floats
    // Somehow the conversion to/from camera coordinates cause precision loss. 
    // This causes the rendered image has grid artifacts in light tracing. 
    // Workaround: use 64-bit float. 
    //
    // However, the following code is slower, so we simply use the above code,
    // which also works well with just 32-bit float.
    /*
	Float x0 = p.x();
	Float y0 = p.y();
	Float z0 = p.z();

	Float xe = eye.x();
	Float ye = eye.y();
	Float ze = eye.z();

    // transform from camera to world
	Vec3 fpZ(0, 0, -focal);	
	fpZ += Vec3(dot(u, eye), dot(v, eye), dot(n, eye));
	fpZ = fpZ.x()*u + fpZ.y()*v + fpZ.z()*n; // correct
    
	Float zf = fpZ.z();	
	Float tem_factor = (zf - ze) / (z0 - ze);
    
	Float _x = tem_factor * (x0 - xe) + xe;
	Float _y = tem_factor * (y0 - ye) + ye;
	Float _z = zf;

	// transform the film point from world space to camera space
	Vec3 f = onb.world_to_local(Vec3(_x, _y, _z));
	f = f - Vec3(dot(u, eye), dot(v, eye), dot(n, eye));
    
	// find the corresponding image pixel
	Vec2 filmPixel(f.x(), f.y());
	pixel = filmPixel * scale + cop;
	*/

    // use truncation to convert continuous coordinates into discrete coordinates (PBRT v2, pg. 338)
    Size2 img_size = sensor->get_image_size();
    int i = (int)(img_size.height - pixel.y());
    int j = (int) pixel.x();
    pixel.set_x(j);
	pixel.set_y(i);

	// check for boundary
	if (i >= 0 && i < img_size.height && j >= 0 && j < img_size.width) {		
		wo_e = unit_vector(eye - p);
		return true;		
	} else {
        // not hit and return a invalid pixel coordinates
		pixel.set_x(-1);
		pixel.set_y(-1);
		return false;
	}	
}

bool Camera::hit(const Vec3 &p, Vec2 &pixel, Vec3 &wo_e, Vec3 &s, Vec3 &sn) {
    // TODO: should discard if p is behind the sensor

	// find the intersection point between the line(eye, p) and the film plane.    
    Vec3 u = onb.u();
    Vec3 v = onb.v();
    Vec3 n = onb.n();    

    // point p in camera coordinates
    Vec3 p0(dot(u, p - eye), dot(v, p - eye), dot(n, p - eye));

    // projected point: P' = E + t (P - E) where P is the 3D point; E is eye.
    // in camera coordinates, E = 0 --> P' = t * P    
    Float t = (-focal) / p0.z();    
    Float xp = t * p0.x();
    Float yp = t * p0.y();
    
    // find the corresponding image pixel
    Vec2 filmPixel2(xp, yp);
    Vec2 cop = sensor->get_cop();
    Vec2 scale = sensor->get_image_film_ratio();
    pixel = filmPixel2 * scale + cop;
    
    // use truncation to convert continuous coordinates into discrete coordinates (PBRT v2, pg. 338)
    Size2 img_size = sensor->get_image_size();
    int i = (int)(img_size.height - pixel.y());
    int j = (int) pixel.x();
    pixel.set_x(j);
	pixel.set_y(i);

	// check for boundary
	if (i >= 0 && i < img_size.height && j >= 0 && j < img_size.width) {		
		wo_e = unit_vector(eye - p);

         // transform film pixel from camera space into world space
        Vec3 filmPixel3(xp, yp, -focal);

        s  = eye + onb.local_to_world(filmPixel3);
        sn = -onb.n();

		return true;		
	} else {
        // not hit and return a invalid pixel coordinates
		pixel.set_x(-1);
		pixel.set_y(-1);
		return false;
	}	
}

/**
 * This is a modification of the above hit() method.
 */
bool Camera::get_sensor_point(const Vec3 &p, Vec3 &s, Vec3 &sn) {
    Vec2 pixel;
    
    Vec3 u = onb.u();
    Vec3 v = onb.v();
    Vec3 n = onb.n();    
        
    Vec3 p0(dot(u, p - eye), dot(v, p - eye), dot(n, p - eye));
        
    Float t = (-focal) / p0.z();
    Float xp = t * p0.x();
    Float yp = t * p0.y();
        
    Vec2 filmPixel2(xp, yp);
    Vec2 cop = sensor->get_cop();
    Vec2 scale = sensor->get_image_film_ratio();
    pixel = filmPixel2 * scale + cop;

    Size2 img_size = sensor->get_image_size();
    int i = (int)(img_size.height - pixel.y());
    int j = (int) pixel.x();
    pixel.set_x(j);
	pixel.set_y(i);

	if (i >= 0 && i < img_size.height && j >= 0 && j < img_size.width) {		
        // transform film pixel from camera space into world space
        Vec3 filmPixel3(xp, yp, -focal);

        s  = eye + onb.local_to_world(filmPixel3);
        sn = -onb.n();

		return true;		
	}  
	return false;		
}

Vec3 Camera::get_sensor_point(const Vec2 &pixel) {
    Vec2 cop = sensor->get_cop();
    Vec2 scale = sensor->get_image_film_ratio();
    Vec2 film_pixel = (pixel - cop) / scale;

    return Vec3(film_pixel.x(), film_pixel.y(), -focal);
}

void Camera::set_default_pixel_list() {        
    pixels.clear();
    Size2 img_size = sensor->get_image_size();
    for (int j = 0; j < img_size.width; ++j) {
        for (int i = 0; i < img_size.height; ++i) {            
            pixels.push_back(Vec2(j, i));
        }
    }
}

void Camera::read_pixel_list(const string file) {    
    FILE *f = fopen(file.c_str(), "r");
    if (f == NULL) return;

    int n, i, j;
    fscanf(f, "%d", &n);
    if (n <= 0) {
        fclose(f);
        return;
    }
    
    pixels.clear();
    for (int k = 0; k < n; ++k) {
        fscanf(f, "%d %d", &i, &j); // row column format
        pixels.push_back(Vec2(j, i)); // j --> x is column
    }
    fclose(f);
}

vector<Vec2> Camera::get_pixel_list() const {
    return pixels;
}

void Camera::set_log_pixel(bool log) {
    is_log_pixel = log;
}

bool Camera::get_log_pixel() const {
    return is_log_pixel;
}

ofstream& Camera::get_log_stream() {
    return log_stream;
}

void Camera::open_log_stream(const string file) {
    log_stream.open(file.c_str());
}

/*
bool Camera::has_pixels() const {
    return sensor->has_pixels();
}

void Camera::next_pixel(int &i, int &j) {    
    pixel_generator->next_pixel(i, j);
}*/

} // end namespace Renzoku
