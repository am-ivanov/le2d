CFLAGS = -fopenmp -O2 -DUSE_SOA 

all: le_core.h le_core.c main.c
	icc $(CFLAGS) -DAUTOVECT -o cf_av main.c le_core.c -lm
	icc $(CFLAGS) -o cf main.c le_core.c -lm
	gcc $(CFLAGS) -o cf_gcc main.c le_core.c -lm
