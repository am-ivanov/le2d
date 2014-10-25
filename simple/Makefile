CFLAGS = -fopenmp -O2

all: le_core.h le_core.c main.c
	gcc $(CFLAGS) -o simple_gcc main.c le_core.c -lm
	icc $(CFLAGS) -o simple main.c le_core.c -lm