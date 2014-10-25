CFLAGS = -fopenmp -O2 -DUSE_SOA

all: le_core.h le_core.c main.c
	gcc $(CFLAGS) -o spl_gr_gcc main.c le_core.c -lm
	icc $(CFLAGS) -o spl_gr main.c le_core.c -lm
	icc $(CFLAGS) -DAUTOVECT -o spl_gr_av main.c le_core.c -lm
