CFLAGS = -pthread -DUSE_SOA -O2

all: le_core.h le_core.c main.c 
	gcc $(CFLAGS) -o pthread_gcc main.c le_core.c -lm
	icc $(CFLAGS) -o pthread main.c le_core.c -lm
	icc $(CFLAGS) -fopenmp -DAUTOVECT -o pthread_av main.c le_core.c -lm
