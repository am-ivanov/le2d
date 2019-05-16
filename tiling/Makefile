CFLAGS = -O2 -g #-fopenmp 

all: gcc_compiler #icc_compiler

gcc_compiler: le_core.h le_core.c main.c
	gcc $(CFLAGS) -o tiled_gcc_f main.c le_core.c -lm
	gcc $(CFLAGS) -DUSE_DOUBLE -o tiled_gcc_d main.c le_core.c -lm
	
#gcc $(CFLAGS) -fno-tree-vectorize -o simple_gcc_f main.c le_core.c -lm
#gcc $(CFLAGS) -fno-tree-vectorize -DUSE_DOUBLE -o simple_gcc_d main.c le_core.c -lm

icc_compiler: le_core.h le_core.c main.c
	icc $(CFLAGS) -fno-tree-vectorize -o simple_icc_f main.c le_core.c -lm
	icc $(CFLAGS) -fno-tree-vectorize -DUSE_DOUBLE -o simple_icc_d main.c le_core.c -lm

clean: 
	rm simple_*
