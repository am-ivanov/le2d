CFLAGS = -pthread -O2

all: gcc_compiler icc_compiler

gcc_compiler: le_core.h le_core.c main.c
	gcc $(CFLAGS) -fno-tree-vectorize -o pthread_gcc_f main.c le_core.c -lm
	gcc $(CFLAGS) -fno-tree-vectorize -DUSE_DOUBLE -o pthread_gcc_d main.c le_core.c -lm

icc_compiler: le_core.h le_core.c main.c
	icc $(CFLAGS) -fno-tree-vectorize -o pthread_icc_f main.c le_core.c -lm
	icc $(CFLAGS) -fno-tree-vectorize -DUSE_DOUBLE -o pthread_icc_d main.c le_core.c -lm
	icc $(CFLAGS) -fopenmp -DAUTOVECT -msse -o pthread_icc_sse_f main.c le_core.c -lm
	icc $(CFLAGS) -fopenmp -DAUTOVECT -DUSE_DOUBLE -msse -o pthread_icc_sse_d main.c le_core.c -lm
	icc $(CFLAGS) -fopenmp -DAUTOVECT -mavx -o pthread_icc_avx_f main.c le_core.c -lm
	icc $(CFLAGS) -fopenmp -DAUTOVECT -DUSE_DOUBLE -mavx -o pthread_icc_avx_d main.c le_core.c -lm

clean: 
	rm pthread_*