#include "base_object.h"
#include <xmmintrin.h>

using namespace Renzoku;

BaseObject::BaseObject(void) {
}


BaseObject::~BaseObject(void) {
}

void *BaseObject::operator new(size_t size) {
    return _mm_malloc(size, 16);
}

void BaseObject::operator delete(void *ptr) {
    _mm_free(ptr);
}

void *BaseObject::operator new[](size_t size) {
    return _mm_malloc(size, 16);
}

void BaseObject::operator delete[](void *ptr) {
    _mm_free(ptr);
}

