CC = nvcc
CFLAGS = 

all: le_core.cu main.cu
	$(CC) $(CFLAGS) le_core.cu main.cu -O2 -o cuda_le2d -lm
	$(CC) $(CFLAGS) le_core.cu main.cu -g -DDEBUG -o cuda_le2d_debug -lm

clean: 
	rm cuda_le2d
