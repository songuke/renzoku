/* predefined class names */
#ifndef _COMMON_H_
#define _COMMON_H_

#define _USE_MATH_DEFINES   // for M_PI to be included with <cmath> in Visual Studio
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>			// min/max
using namespace std;

#ifdef _MSC_VER 
	#define OPENEXR_DLL         // for half::_eLut and half::_toFloat linking http://hebbut.net/Public.Offerings/OpenEXR.html
	#define NOMINMAX            // disable min/max in windefs.h
#endif

#define __align16 __declspec(align(16)) // For alignment of 16 bytes during static memory allocation in stack.

#define FREEIMAGE

/* Ensure automatic linking by Boost looks for 'boost-...-.lib' instead of 'libboost-...-.lib'. 
 *
 * Current bug: 
 * Thread will look for libboost-... while unit test looks for boost-...
 * They are all dynamic libs but the naming to look for are not consistent.
 */
#define BOOST_ALL_DYN_LINK

namespace Renzoku {
//#define _MATH_SSE_              // comment out this to revert to old vector library

typedef int             Int;
typedef float           Float;            
typedef unsigned char   Byte;

typedef int             ID;

const Float A_PI    = 3.14159265358979f;
const Float INV_PI  = 0.31830988618379f;  // 14 digits
const Float TWO_PI  = 6.28318530717959f;
const Float HALF_PI = 1.57079632679490f;
const Float INV_2PI = 0.15915494309190f;
const Float INV_4PI = 0.07957747154595f;

const Float RADIAN_TO_DEGREE = 180.0f / A_PI;
const Float DEGREE_TO_RADIAN = A_PI / 180.0f;

const Float RAY_TMIN = 1e-3f;             // Setting this to 1e-5 can cause "self-intersection" (PBRT 3.1.4). 
                                          // Test with single uniform eye ray and direct lighting (fixed point light or directional light) 
const Float RAY_TMAX = 1e6f;
const int   RAY_BUNDLE_DIM  = 64;
const int   RAY_BUNDLE_SIZE = RAY_BUNDLE_DIM * RAY_BUNDLE_DIM;       // DIM^2 = SIZE

const Float ONB_EPSILON     = 0.01f;
const Float ZERO_EPSILON    = 1e-6f;

template <class T> class Image;

typedef Image<Byte>		ImageByte;
typedef Image<Int>      ImageInt;
typedef Image<Float>	ImageFloat;

class BaseObject;
class Vec2;
class Vec3;
class Rgb;
class Onb;
class Ray;
class Random;

struct Size2;
struct Size3;

class HitRecord;
struct Receiver;
struct LocalGeometry;

class Shape;
class Bsdf;
class Material;
class Texture;
class Surface;
class Aggregate;
class Scene;

class Frame;
class Renderer;
class Sampler;

class Integrator;
class FrameBuffer;      // to replace Frame

class Camera;
class CameraView;

class Scheduler;
class ThreadController;

class Light;
class AreaLight;
class PointLight;
class SpotLight;
class BrdfPointLight;
class EnvLight;
class DirectionalLight;

class Triangle;
class Sphere;

class Stats;

class Viewer;
class MeshView;
class ImageBlockView;

typedef vector<Shape*>              Shapes;                     // store a list of pointers for polymorphism
typedef vector<Material*>           Materials;                  // to achieve contiguous memory, use custom allocator, or
                                                                // store many arrays, each for an object type (discard polymorphism)
typedef vector<Texture*>            Textures;

typedef vector<Light*>              Lights;                     // BVH stores a local contiguous copy of triangles for faster intersection
typedef vector<BrdfPointLight*>     BrdfPointLightPtrs;
typedef vector<BrdfPointLight>      BrdfPointLights;

typedef vector<Camera *>            Cameras;
typedef vector<CameraView *>        CameraViews;

};

#endif
