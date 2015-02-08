/*
 * Author: Nikolay Khokhlov <k_h@inbox.ru>, (C) 2012
 * Compile: gcc -O2 -o le2d *.c -lm
 */

#include <stdio.h>
#include <mpi.h>

#include "le_core.h"

#ifdef SAVE_EVERY_STEPS
const int save_every_step = 1;
#else
const int save_every_step = 0;
#endif

int main(int argc, char *argv[])
{
	int i;
	le_point2 n;
	int steps;
	le_task task;
	le_material mat;
	le_vec2 h;
	real dt;
	le_vec2 center;
	
	char name[1000];
	int size, rank;
	
	double t;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (4 != argc) {
		if (0 == rank)
			printf("Usage: %s nx ny steps.\n", argv[0]);
		MPI_Finalize();
		return 0;
	}
	
	n.x = atoi(argv[1]);
	n.y = atoi(argv[2]);
	steps = atoi(argv[3]);
	h.x = 1.0; 
	h.y = 1.0;
	dt = 0.3;
	center.x = n.x / 2;
	center.y = n.y / 2;

	le_init_material(2.0, 1.0, 1.5, &mat);
	le_init_task(&task, dt, h, mat, n);
	if (0 == rank) {
		le_init_grid(&task);
		le_set_ball(&task, center, 10.0, 1.0);
	}
	scatter_work(&task, size, rank);
	//if (0 == rank) printf("ok\n");
	MPI_Barrier(MPI_COMM_WORLD);
	t = MPI_Wtime();
	for (i = 0; i < steps; i++) {
		// WARNING: le_save_task is not working here
		if (save_every_step) {
			sprintf(name, "out-%06d.vtk", i);
			le_save_task(&task, name);
		}
		le_step(&task);	
	}
	MPI_Barrier(MPI_COMM_WORLD);
	t = MPI_Wtime() - t;
	
	if (0 == rank)
		printf("%d %d %d %f\n", n.x, n.y, steps, t);
	gather_work(&task, size);
	if (0 == rank)
		le_save_task(&task, "result.vtk");
	if (0 == rank)
		le_free_task(&task);
	MPI_Finalize();
	return 0;
}
