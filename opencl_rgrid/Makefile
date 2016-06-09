CC = g++
CFLAGS = -O2

OPENCL_INC = /usr/local/cuda/include
OPENCL_LIB = /usr/local/cuda/lib64

CUR_DIR = ./

#OPENCL_INC = /opt/AMDAPPSDK-3.0-0-Beta/include
#OPENCL_LIB = /opt/AMDAPPSDK-3.0-0-Beta/lib/x86_64/

all: le_core.c main.c
	$(CC) $(CFLAGS) -I $(CUR_DIR) -L $(CUR_DIR) main.c -c -o main_f.o
	$(CC) $(CFLAGS) -I $(CUR_DIR) -L $(CUR_DIR) -I $(OPENCL_INC) le_core.c -c -o le_core_f.o
	$(CC) $(CFLAGS) -I $(CUR_DIR) -L $(CUR_DIR) -L $(OPENCL_LIB) main_f.o le_core_f.o -o le2d_opencl_f -lOpenCL -lrgrid -lm
	$(CC) $(CFLAGS) -I $(CUR_DIR) -L $(CUR_DIR) main.c -DUSE_DOUBLE -c -o main_d.o
	$(CC) $(CFLAGS) -I $(CUR_DIR) -L $(CUR_DIR) -I $(OPENCL_INC) le_core.c -DUSE_DOUBLE -c -o le_core_d.o
	$(CC) $(CFLAGS) -I $(CUR_DIR) -L $(CUR_DIR) -L $(OPENCL_LIB) main_d.o le_core_d.o -DUSE_DOUBLE -o le2d_opencl_d -lOpenCL -lrgrid -lm

clean: 
	rm *.o le2d_opencl_*
