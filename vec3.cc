#include "common.h"

#ifdef _MATH_SSE_
#include "vec3_sse.cc"
#else
#include "vec3_cpu.cc"
#endif
