#ifndef _CAMERA_H_
#define _CAMERA_H_

#include "common.h"
#include "named_object.h"
#include "base_object.h"
#include "math3.h"
#include "ray.h"
#include "sensor.h"
#include "boundingbox.h"

#include <fstream>
using namespace std;

namespace Renzoku {


/**
 * Pin-hole camera model
 *   Right-handed, look at negative Z-axis.
 *   Image (pixel) plane origin: top left.
 *   Film plane origin: bottom left (following right-handed coordinates).
 *
 * World to camera matrix:
 * [ u * * -dot(u, eye)
 *   v * * -dot(v, eye)
 *   n * * -dot(n, eye)
 *   0 0 0 1            ]
 *
 * Camera to world matrix:
 * [ u v n  dot(u, eye) 
 *   * * *  dot(v, eye) 
 *   * * *  dot(n, eye)
 *   0 0 0  1           ]
 *
 * 
 * Camera fetches a SensorSample and turns it into a SensorRay.
 * Therefore, it maps a SensorSamplePack into a RayBundle.
 * Scheduler is responsible for scheduling the processor that calculates information (RGB, indirect bounces, depth, normal, etc.)
 * that each ray bundle needs.
 * 
 */
class Camera : public BaseObject, public NamedObject {
public:
    Camera();

    Camera(const Vec3 &_eye, const Vec3 &_lookat, const Vec3 &_up, Float _focal, Sensor *sensor);
    
    Camera *clone();

    /**
     * Set camera parameters based on settings from OpenGL camera.
     */
    void set_perspective(const Vec3 &_eye, const Vec3 &_lookat, const Vec3 &_up, Float fov);

    void set_perspective(const Vec3 &_eye, const Vec3 &_lookat, const Vec3 &_up, Float _focal, 
                         Sensor *sensor);
    // pass Vec3 in by reference to avoid 16-byte alignment problem when using the SSE __m128 struct in Vec3.

    inline Size2 get_film_size() const;
    inline Size2 get_image_size() const;
    
    inline Sensor *get_sensor() const;
    inline void set_sensor(Sensor *sensor);

    /**
     * Update center of projection every time sensor film/image size changes.
     */
    inline void boot_sensor();

    inline Onb get_onb() const;

    /**
     * Return vertical field of view of the camera. For use with gluPerspective in OpenGL. 
     *
     * Fov is returned in degree. 
     */
    inline Float get_vertical_fov() const;

    /**
     * Set the focal length of the lens based on vertical field of view and the height of the film size. 
     *
     * Field of view is defined in degree unit. 
     */
    inline void set_vertical_fov(Float fov); 
    
    inline Float get_focal_length() const;

    /**
     * Utility methods to convert between vertical field of view (OpenGL) and focal length (physical lens)
     */
    static Float focal_length_from_fov(Float fov, Float film_size_height = 0.025f);
    static Float fov_from_focal_length(Float focal, Float film_size_height = 0.025f);

    inline Vec3 get_eye()    const;
    inline Vec3 get_lookat() const;
    inline Vec3 get_view_dir() const;
    inline Vec3 get_up()     const;
    inline Float get_near_plane() const;
    inline Float get_far_plane() const;
    inline void set_near_plane(Float near_plane);
    inline void set_far_plane(Float far_plane);

    /**
     * Generate a ray from the camera origin through a pixel on the image plane.
     */
    Ray shoot_ray(int i, int j, Random &rd);
    Ray shoot_ray(int i, int j, Vec2 jitter);

    /**
     * Generate a ray from the camera origin through a pixel center (jitter 0.5) on the image plane.
     */
    Ray shoot_ray(int i, int j);
    Ray shoot_ray(const SensorSample &sample);

	/*
	 * Check if a ray from the scene hits the sensor.
	 * If yes, return the pixel coordinates in a Vec2 and the direction from 
     * p to eye
     */
	bool hit(const Vec3 &p, Vec2 &pixel, Vec3 &wo_e);	

    /**
     * Check if a ray from the scene hits the sensor.
	 * If yes, return the pixel coordinates in a Vec2 and the direction from 
     * p to eye, and the hit point and its normal on the sensor.
     */
    bool hit(const Vec3 &p, Vec2 &pixel, Vec3 &wo_e, Vec3 &s, Vec3 &sn);

    /**
     * Check if a ray from the scene hits the sensor.
	 * If yes, return the hit point and its normal on the sensor 
     * 
     * p: the visible surface point.
     */
    bool get_sensor_point(const Vec3 &p, Vec3 &s, Vec3 &sn);

    /**
     * Return the corresponding point on the film plane (in camera coordinate frame)
     * Return value: [x, y, z] where z is -focal.
     */
    Vec3 get_sensor_point(const Vec2 &pixel);

    void read_pixel_list(const string file);
    
    vector<Vec2> get_pixel_list() const;
    void set_default_pixel_list();
    void set_log_pixel(bool log);
    bool get_log_pixel() const;

    void open_log_stream(const string file);
    ofstream& get_log_stream();

    /**
     * Return the probability of the visible surface point obtained 
     * by intersecting a ray from eye through a point sampled in a pixel.
     */
    Float pdf(const Vec3 &p, const Vec3 &n);

public:
    /**
     * Approximate clipping planes for shadow map projection
     */
    static void approximate_near_far(const Vec3 &p, const Vec3 &n, const BoundingBox &bb, 
                                     Float &nn, Float &ff) {
    
        Onb uvn;
        uvn.init_from_n(n);
    
        Vec3 ma = uvn.world_to_local(bb.v_min - p);
        Vec3 mb = uvn.world_to_local(Vec3(bb.v_min.x(), bb.v_min.y(), bb.v_max.z()) - p);
        Vec3 mc = uvn.world_to_local(Vec3(bb.v_min.x(), bb.v_max.y(), bb.v_min.z()) - p);
        Vec3 md = uvn.world_to_local(Vec3(bb.v_min.x(), bb.v_max.y(), bb.v_max.z()) - p);
        Vec3 me = uvn.world_to_local(Vec3(bb.v_max.x(), bb.v_min.y(), bb.v_min.z()) - p);
        Vec3 mf = uvn.world_to_local(Vec3(bb.v_max.x(), bb.v_min.y(), bb.v_max.z()) - p);
        Vec3 mg = uvn.world_to_local(Vec3(bb.v_max.x(), bb.v_max.y(), bb.v_min.z()) - p);
        Vec3 mh = uvn.world_to_local(bb.v_max - p);
   
        // the view frustum is along n (for a light we look at positive z-axis)
        nn = std::max(0.1f, -1.0f + std::min(
                                std::min(std::min(ma.z(), mb.z()), std::min(mc.z(), md.z())),
                                std::min(std::min(me.z(), mf.z()), std::min(mg.z(), mh.z()))));
        
        ff = std::max(0.1f, 1.0f + std::max(
                                std::max(std::max(ma.z(), mb.z()), std::max(mc.z(), md.z())),
                                std::max(std::max(me.z(), mf.z()), std::max(mg.z(), mh.z()))));                            
    }

protected:
    Vec3 eye, lookat, up;       // extrinsic parameters
    Float focal;    
    
    Onb onb;

    Float near_plane, far_plane; // for rasterization use (e.g., OpenGL, VPL on GPU)

    Sensor *sensor;    
        
    vector<Vec2> pixels;        // a hack for collecting pixel values for inspection
    bool is_log_pixel;
    ofstream log_stream;

    Scene *scene;
};

inline Size2 Camera::get_image_size() const { 
    return sensor->get_image_size();
}

inline Size2 Camera::get_film_size() const { 
    return sensor->get_film_size(); 
}

inline Onb Camera::get_onb() const {
    return onb;
}

inline Float Camera::get_vertical_fov() const {   
    Size2 film_size = sensor->get_film_size();
    return 2 * atan(film_size.height / (2.0f * focal)) * RADIAN_TO_DEGREE;
}

inline void Camera::set_vertical_fov(Float fov) {
    Size2 film_size = sensor->get_film_size();
    focal = film_size.height / (2.0f * tan(fov / 2 * DEGREE_TO_RADIAN));
}

inline Float Camera::focal_length_from_fov(Float fov, Float film_size_height) {
    return film_size_height / (2.0f * tan(fov / 2 * DEGREE_TO_RADIAN));
}

inline Float Camera::fov_from_focal_length(Float focal, Float film_size_height) {
    return 2 * atan(film_size_height / (2.0f * focal)) * RADIAN_TO_DEGREE;
}

inline Float Camera::get_focal_length() const {
    return focal;
}

inline Vec3 Camera::get_eye() const {
    return eye;
}

inline Vec3 Camera::get_lookat() const {
    return lookat;
}

inline Vec3 Camera::get_view_dir() const {
    return unit_vector(lookat - eye);
}

inline Vec3 Camera::get_up() const {
    return up;
}

inline void Camera::set_near_plane(Float near_plane) {
    this->near_plane = near_plane;
}

inline Float Camera::get_near_plane() const {
    return near_plane;
}

inline void Camera::set_far_plane(Float far_plane) {
    this->far_plane = far_plane;
}

inline Float Camera::get_far_plane() const {
    return far_plane;
}

inline Sensor* Camera::get_sensor() const {
    return sensor;
}

inline void Camera::set_sensor(Sensor *sensor) {
    this->sensor = sensor;
}

} // end namespace Renzoku 

#endif

