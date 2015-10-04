CC = /usr/local/cuda-6.5/bin/nvcc

all: le_core.cu main.cu le_core.h
	$(CC) $(CFLAGS) le_core.cu main.cu -DUSE_DOUBLE -O2 -o cuda_le2d_d -lm
	$(CC) $(CFLAGS) le_core.cu main.cu -O2 -o cuda_le2d_f -lm

debug: le_core.cu main.cu le_core.h
	$(CC) $(CFLAGS) le_core.cu main.cu -g -DDEBUG -o cuda_le2d_debug_f -lm
	$(CC) $(CFLAGS) le_core.cu main.cu -DUSE_DOUBLE -g -DDEBUG -o cuda_le2d_debug_d -lm
clean: 
	rm cuda_le2d*
