// This file include OpenGL headers without GLUT. 
// Code ripped from freeglut_std.h.

#ifdef __cplusplus
    extern "C" {
#endif

#ifdef _WIN32

    #  ifndef WIN32_LEAN_AND_MEAN
    #    define WIN32_LEAN_AND_MEAN 1
    #  endif
    // min and max conflicts will cause several non-sense error http://support.microsoft.com/kb/143208
    #  ifndef NOMINMAX
    #    define NOMINMAX
    #  endif 
    // windows.h is needed for WINGDIAPI keyword
    #include <windows.h>
#else
    #define APIENTRY
#endif

#if __APPLE__
#   include <OpenGL/gl.h>
#   include <OpenGL/glu.h>
#else
#   include <GL/gl.h>
#   include <GL/glu.h>
#endif

#ifdef __cplusplus
    }
#endif
