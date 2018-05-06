#ifndef _SIZE_H_
#define _SIZE_H_

#include "vec2.h"
#include "vec3.h"

namespace Renzoku {

struct Size2 {
	Float width;
	Float height;
		
    Size2() : width(0), height(0) {}
	Size2(Float width, Float height) : width(width), height(height) {}
	inline Size2 operator*(Float v) const { 
		return Size2(width * v, height * v);
	}
	
	inline Size2 operator*(Vec2 v) const {
		return Size2(width * v.x(), height * v.y());
	}
	
	inline Size2 operator*(Size2 v) const {
		return Size2(width * v.width, height * v.height);
	}
	
	inline Size2 operator/(Float v) const { 
		Float inv_v = 1.0f / v;
		return Size2(width * inv_v, height * inv_v);
	}
	
	inline Size2 operator/(Vec2 v) const { 		
		return Size2(width / v.x(), height / v.y());
	}
	
	inline Size2 operator/(Size2 v) const { 		
		return Size2(width / v.width, height / v.height);
	}

    inline bool operator==(Size2 v) const {
        if ((height == v.height) && (width == v.width)) 
            return true;
        return false;
    }

    inline Vec2 to_vec2() const {
        return Vec2(width, height);
    }
};

struct Size3 {
	Float width;
	Float height;
	Float depth;
	
    Size3() : width(0), height(0), depth(0) {}
	Size3(Float width, Float height, Float depth) : width(width), height(height), depth(depth) {}
	inline Size3 operator*(Float v) const { 
		return Size3(width * v, height * v, depth * v);
	}
	
	inline Size3 operator/(Float v) const { 
		Float inv_v = 1.0f / v;
		return Size3(width * inv_v, height * inv_v, depth * inv_v);
	}

    inline Vec3 to_vec3() const {
        return Vec3(width, height, depth);
    }
};

} // end namespace

#endif