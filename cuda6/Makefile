CC = nvcc

all: le_core.cpp main.cpp le_core.h
	$(CC) $(CFLAGS) le_core.cpp main.cpp -x cu -DUSE_DOUBLE -O2 -o cuda_le2d_d -lm
	$(CC) $(CFLAGS) le_core.cpp main.cpp -x cu -O2 -o cuda_le2d_f -lm

debug: le_core.cpp main.cpp le_core.h
	$(CC) $(CFLAGS) le_core.cpp main.cpp -g -DDEBUG -o cuda_le2d_debug_f -lm
	$(CC) $(CFLAGS) le_core.cpp main.cpp -DUSE_DOUBLE -g -DDEBUG -o cuda_le2d_debug_d -lm
clean: 
	rm cuda_le2d*
