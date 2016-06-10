OPENCL_INC = /usr/local/cuda/include
OPENCL_LIB = /usr/local/cuda/lib64

CXXFLAGS = -O2 -I. -L. -I $(OPENCL_INC) -L $(OPENCL_LIB) -DUSE_OPENCL

#OPENCL_INC = /opt/AMDAPPSDK-3.0-0-Beta/include
#OPENCL_LIB = /opt/AMDAPPSDK-3.0-0-Beta/lib/x86_64/

all: le2d_opencl_f le2d_opencl_d

le2d_opencl_f: $(RGRID_OBJ_FILES)
	$(CXX) $(CXXFLAGS) $(wildcard rgrid/src/*.cpp) le_core.c main.c -o le2d_opencl_f -lOpenCL -lm

le2d_opencl_d: $(RGRID_OBJ_FILES)
	$(CXX) $(CXXFLAGS) $(wildcard rgrid/src/*.cpp) le_core.c main.c -DUSE_DOUBLE -o le2d_opencl_d -lOpenCL -lm

.PHONY: clean 
clean:
	rm -f *.o le2d_opencl_*
