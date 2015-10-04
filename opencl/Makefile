CC = gcc
CFLAGS = -O2

OPENCL_INC = /usr/local/cuda/include
OPENCL_LIB = /usr/local/cuda/lib64

#OPENCL_INC = /opt/AMDAPPSDK-3.0-0-Beta/include
#OPENCL_LIB = /opt/AMDAPPSDK-3.0-0-Beta/lib/x86_64/

all: le_core.c main.c
	$(CC) $(CFLAGS) main.c -c -o main_f.o
	$(CC) $(CFLAGS) -I $(OPENCL_INC) le_core.c -c -o le_core_f.o
	$(CC) $(CFLAGS) -L $(OPENCL_LIB) main_f.o le_core_f.o -o le2d_opencl_f -lOpenCL -lm
	$(CC) $(CFLAGS) main.c -DUSE_DOUBLE -c -o main_d.o
	$(CC) $(CFLAGS) -I $(OPENCL_INC) le_core.c -DUSE_DOUBLE -c -o le_core_d.o
	$(CC) $(CFLAGS) -L $(OPENCL_LIB) main_d.o le_core_d.o -DUSE_DOUBLE -o le2d_opencl_d -lOpenCL -lm

clean: 
	rm *.o le2d_opencl_*
