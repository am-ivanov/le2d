CFLAGS = -O2 -std=c++14#-fopenmp

all: gcc_compiler #icc_compiler

gcc_compiler: le_core.h le_core.cpp main.cpp gen.h
	g++ $(CFLAGS) -o tiled_gcc_f main.cpp le_core.cpp -lm
	g++ $(CFLAGS) -DUSE_DOUBLE -o tiled_gcc_d main.cpp le_core.cpp -lm

#gcc $(CFLAGS) -fno-tree-vectorize -o simple_gcc_f main.c le_core.c -lm
#gcc $(CFLAGS) -fno-tree-vectorize -DUSE_DOUBLE -o simple_gcc_d main.c le_core.c -lm

#icc_compiler: le_core.h le_core.cpp main.cpp
#	icc $(CFLAGS) -fno-tree-vectorize -o simple_icc_f main.c le_core.c -lm
#	icc $(CFLAGS) -fno-tree-vectorize -DUSE_DOUBLE -o simple_icc_d main.c le_core.c -lm

clean:
	rm simple_*
