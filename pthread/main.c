/*
 * Author: Nikolay Khokhlov <k_h@inbox.ru>, (C) 2012
 * Compile: gcc -O2 -o le2d *.c -lm
 */

#include <stdio.h>
#include <pthread.h>
#include <assert.h>

#include "le_core.h"

#ifdef SAVE_EVERY_STEPS
const int save_every_step = 1;
#else
const int save_every_step = 0;
#endif

#define NUM_THREADS 100

extern pthread_barrier_t b;

int main(int argc, char *argv[])
{
	if (argc != 5) {
		printf("Usage: %s nx ny steps threads.\n", argv[0]);
		return 1;
	}
	int i, ti;
	le_point2 n = {atoi(argv[1]), atoi(argv[2])};
	int steps = atoi(argv[3]);
	int max_threads = atoi(argv[4]);
	le_task task;
	le_material mat;
	le_vec2 h = {1.0, 1.0};
	real dt = 0.3;
	le_vec2 center = {n.x / 2, n.y / 2};
	char name[1000];

	double t;
	unsigned long cc;
	
	/* pthread variables */
	pthread_t threads[NUM_THREADS];
	st_pthread data[NUM_THREADS];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	/* end of pthread variables */

	le_init_material(2.0, 1.0, 1.5, &mat);
	task.max_threads = max_threads;
	task.steps = steps;
	le_init_task(&task, dt, h, mat, n);
	le_set_ball(&task, center, 10.0, 1.0);

	pthread_barrier_init(&b, NULL, max_threads);
	
	cc = getCC();
	t = timer();
	
	for(ti = 0; ti < max_threads; ti++) {		
		data[ti].task = &task;
		data[ti].thread_num = ti;
		int rc = pthread_create(&(threads[ti]), &attr, le_step, (void *)(&(data[ti])));
		assert(!rc);
	}
	for(ti = 0; ti < max_threads; ti++) {
		void *st;
		int rc = pthread_join(threads[ti], &st);
		assert(!rc);;
	}
	
	t = timer() - t;
	cc = getCC() - cc;
	
	pthread_barrier_destroy(&b);
	
	printf("%d %d %d %f ", n.x, n.y, steps, t);
	
	le_save_task(&task, "result.vtk");

	le_free_task(&task);
	pthread_exit(NULL);
}
