#ifndef _MAT4_H_
#define _MAT4_H_

#include "common.h"

namespace Renzoku {

class Mat4 {
public:
	Mat4();
	Mat4(const Mat4 &a);
	Mat4(const Float m[16]);
	
	Float* get_data();
	
	Float  operator[](int i) const;
	Float& operator[](int i);
	Float  operator()(int i, int j) const;
	Float& operator()(int i, int j);
	
    Mat4& operator=(const Mat4 &a);
    Mat4 operator*(const Mat4 &a) const;

protected:
	Float m[16];		// column major
};

inline Mat4::Mat4() {
	m[0] = 0;
	m[1] = 0;
	m[2] = 0;
	m[3] = 0;
	m[4] = 0;
	m[5] = 0;
	m[6] = 0;
	m[7] = 0;
	m[8] = 0;
	m[9] = 0;
	m[10] = 0;
	m[11] = 0;
	m[12] = 0;
	m[13] = 0;
	m[14] = 0;
    m[15] = 0;
}

inline Mat4::Mat4(const Float a[16]) {
	m[0] = a[0];
	m[1] = a[1];
	m[2] = a[2];
	m[3] = a[3];
	m[4] = a[4];
	m[5] = a[5];
	m[6] = a[6];
	m[7] = a[7];
	m[8] = a[8]; 
	m[9] = a[9]; 
	m[10] = a[10];
	m[11] = a[11];
	m[12] = a[12];
	m[13] = a[13];
	m[14] = a[14];
	m[15] = a[15];
}

inline Mat4::Mat4(const Mat4 &a) {
	m[0] = a.m[0];
	m[1] = a.m[1];
	m[2] = a.m[2];
	m[3] = a.m[3];
	m[4] = a.m[4];
	m[5] = a.m[5];
	m[6] = a.m[6];
	m[7] = a.m[7];
	m[8] = a.m[8]; 
	m[9] = a.m[9]; 
	m[10] = a.m[10];
	m[11] = a.m[11];
	m[12] = a.m[12];
	m[13] = a.m[13];
	m[14] = a.m[14];
	m[15] = a.m[15];
}

inline Float* Mat4::get_data() {
	return &m[0];
}

inline Float Mat4::operator[](int i) const {
	return m[i];
}

inline Float& Mat4::operator[](int i) {
	return m[i];
}

inline Float Mat4::operator()(int i, int j) const {
	return m[j * 4 + i];
}
	
inline Float& Mat4::operator()(int i, int j) {
	return m[j * 4 + i];
}

inline Mat4 Mat4::operator*(const Mat4 &a) const {
    Mat4 out;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 4; ++k) {
                out(i, j) += (*this)(i, k) * a(k, j);
            }
        }
    }
    return out;
}

inline Mat4& Mat4::operator=(const Mat4 &a) {
    m[0] = a.m[0];
	m[1] = a.m[1];
	m[2] = a.m[2];
	m[3] = a.m[3];
	m[4] = a.m[4];
	m[5] = a.m[5];
	m[6] = a.m[6];
	m[7] = a.m[7];
	m[8] = a.m[8]; 
	m[9] = a.m[9]; 
	m[10] = a.m[10];
	m[11] = a.m[11];
	m[12] = a.m[12];
	m[13] = a.m[13];
	m[14] = a.m[14];
	m[15] = a.m[15];
    return (*this);
}

} // end namespace

#endif