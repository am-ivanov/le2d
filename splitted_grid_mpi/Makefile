CFLAGS = -O2

all: mpi

mpi: le_core.h le_core.c main.c
	mpicc $(CFLAGS) -DUSE_DOUBLE -o le2d_mpi main.c le_core.c -lm
	
debug:
	mpicc $(CFLAGS) -DUSE_DOUBLE -g -o le2d_mpi main.c le_core.c -lm

clean: 
	rm le2d_mpi