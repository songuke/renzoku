# use Macports' gcc48 since clang on OS X 10.9 fails. The default gcc and g++ still points to clang.
# gcc_select can be used to select gcc if multiple gcc are installed by Macports (does not count OS X's gcc).
GCC = /opt/local/bin/gcc-mp-4.8
GXX = /opt/local/bin/g++-mp-4.8
#GCC = clang
#GXX = clang++
CC = $(GCC)
CPP = $(GXX)

CFLAGS += -Wall -g -O2
# don't enable -ffast-math as it may screw up IEEE 754 standard and the isnan test.

# Linux or Macports
UNAME := $(shell uname -s)
ifeq ($(UNAME), Linux)
USR=/usr
endif
ifeq ($(UNAME), Darwin)
USR=/opt/local
# ensure pkg-config can be found
export PATH := /opt/local/bin:$(PATH)
endif


USR_INC   = -I$(USR)/include/

OPENEXR_INC = -I$(USR)/include/OpenEXR
OPENEXR_LIB = `pkg-config --libs OpenEXR` # -pthread -lIlmImf -lz -lImath -lHalf -lIex -lIlmThread

#MAGICK_INC = -I$(USR)/include/ImageMagick-6
#MAGICK_LIB = `pkg-config --libs Magick++` # -lMagick++ -lMagickCore

ifeq ($(UNAME), Linux)
GLUT_LIB = -lGL -lGLU -lglut
endif
ifeq ($(UNAME), Darwin)
# Qt will link to Apple's OpenGL (AGL) by default. 
# The Qt viewer cannot swap buffers properly if it is linked with OpenGL from Macports (-lGL -lGLU). 
# Use OpenGL framework instead. Use only freeglut from Macports
GLUT_LIB = -framework OpenGL -lglut
endif

GLEW_LIB = `pkg-config --libs glew`

BOOST_LIB = -lboost_system-mt -lboost_thread-mt -lboost_unit_test_framework-mt
FREEIMAGE_LIB = -lFreeImage

EIGEN_INC = -I./vc/3rdparty/eigen

QTDIR = /Users/sonhua/Qt/5.2.1/clang_64
# Qt on OS X is packed as framework
QT_INC = -F$(QTDIR)/lib
QT_LIB = -F$(QTDIR)/lib -framework QtCore -framework QtGui -framework QtWidgets -framework QtOpengl

JSON_INC = -I./vc/3rdparty/libjson
JSON_LIB = -L./vc/3rdparty/libjson -ljson

CUDA_HOME = /Developer/NVIDIA/CUDA-6.0
NVCC = nvcc
CUDA_INC = -I$(CUDA_HOME)/samples/common/inc
CUDA_LIB = -lcudart

COMMON_INC = $(USR_INC) $(MAGICK_INC) $(OPENEXR_INC) $(EIGEN_INC) $(QT_INC) $(JSON_INC)
COMMON_LIB = $(OPENEXR_LIB) $(GLUT_LIB) $(GLEW_LIB) $(FREEIMAGE_LIB) $(BOOST_LIB) $(QT_LIB) $(JSON_LIB)

COMMON_OBJ = error.o vec2.o vec3.o matrix.o matrix_c.o rgb.o base_object.o resource.o log.o stats.o param_dict.o \
             sphere.o triangle.o quad.o surface.o \
             aggregate.o octree.o aggregate_octree.o aggregate_bvh.h aggregate_bvh.cc \
             onb.o sampler.o stratified_sampler.o \
             material.o texture.o diffuse.o thin_transparent.o glass.o mirror.o phong.o ward.o \
             pdf.o \
             point_light.o brdf_point_light.o area_light.o photon_light.o env_light.o dir_light.o spot_light.o \
             sensor.o camera.o \
             mtwist.o \
             image_sampler.o \
             light_particle.o kdtree.o \
             scene.o \
             string_extra.o obj_parser.o obj_list.o obj_loader.o obj_importer.o json_importer.o \
             integrator.o monte_carlo_integrator.o monte_carlo_integrator_impl.o \
             scheduler.o thread_controller.o \
             shading_geometry.o \
             first_bounce.o path_tracing.o vpl.o vpl_segment.o \
             vpl_generator.o bidir_vpl_generator.o vpl_radius_setter.o vpl_density_estimator.o \
             gl_utils.o glsl_shader.o gl_scene.o vpl_gpu.o vsl_gpu.o vpl_segment_gpu.o vpl_adaptive_sampling_gpu.o \
             vpl_ppm.o \
             lightcuts.o cone.o \
             light_transport_matrix.o form_factor_matrix.o brdf_light_matrix.o brdf_surface_matrix.o visibility_matrix.o cluster_builder.o mrcs.o \
             viewer.o glut_viewer.o qt_viewer.o mesh_view.o image_block_view.o gpu_view.o tonemap.o frame.o trackball.o \
             light_tracing.o bidir_path_tracing.o \
             vps.o \
             #cuda/bundle_path_tracing.o cuda/bundle_vpl.o cuda/linear_bvh.o \

INCLUDES = $(wildcard *.h) 
CUDA_INCLUDES = $(wildcard *.cuh)

all: renzoku 

renzoku: renzoku.o $(COMMON_OBJ)
	$(CPP) -o $@ renzoku.o $(COMMON_OBJ) $(COMMON_LIB)

# TODO: auto compile for all tests
#unit_test_vector: mtwist.o randistrs.o stats.o
#	$(CPP) -o $@ unit_test_vector.cc mtwist.o randistrs.o stats.o $(CFLAGS) -I./vc/3rdparty/vectormath/include/sse

test_cone: test_cone.o $(COMMON_OBJ)
	$(CPP) -o $@ test_cone.o $(COMMON_OBJ) $(COMMON_LIB)

test_cone.o: test/test_rzrt/test_cone.cc ${INCLUDES}
	$(CPP) $(CFLAGS) $(COMMON_INC) -I. -c -o $@ $<

%.o: %.cc ${INCLUDES}
	$(CPP) $(CFLAGS) $(COMMON_INC) -c -o $@ $<

%.o: %.cpp ${INCLUDES}
	$(CPP) $(CFLAGS) $(COMMON_INC) -c -o $@ $<

%.o: %.c ${INCLUDES}
	$(CC) $(CFLAGS) $(COMMON_INC) -c -o $@ $<

#%.o: %.cu ${INCLUDES} ${CUDA_INCLUDES}
#	$(NVCC) -ccbin=$(GCC) $(CUDA_INC) -Xcompiler '$(COMMON_INC)' -c -o $@ $<

clean:
	rm *.o

