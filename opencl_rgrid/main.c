/*
 * Author: Nikolay Khokhlov <k_h@inbox.ru>, (C) 2012
 * Compile: gcc -O2 -o le2d *.c -lm
 */

#include <stdio.h>
#include <stdlib.h>

#include "le_types.h"
#include "le_core.h"

/*#ifdef SAVE_EVERY_STEPS
const int save_every_step = 1;
#else
const int save_every_step = 0;
#endif*/

int main(int argc, char *argv[])
{
	if (argc != 4) {
		printf("Usage: %s nx ny steps.\n", argv[0]);
		return 1;
	}
	//int i;
	le_point2 n = {atoi(argv[1]), atoi(argv[2])};
	int steps = atoi(argv[3]);
	le_task task;
	le_material mat;
	le_vec2 h = {1.0, 1.0};
	real dt = 0.3;
	le_vec2 center;
	center.x = n.x / 2;
	center.y = n.y / 2;
	//char name[1000];
	double t;

	le_init_material(2.0, 1.0, 1.5, &mat);
	le_init_task(&task, dt, h, mat, n);
	le_set_ball(&task, center, 10.0, 1.0);

	/*for (i = 0; i < steps; i++) {
		if (save_every_step) {
			sprintf(name, "out-%06d.vtk", i);
			le_save_task(&task, name);
		}
		le_step(&task);
	}*/
	t = le_step(&task, steps);

	printf("%d %d %d %f\n", n.x, n.y, steps, t);

	le_save_task(&task, "result.vtk");

	le_free_task(&task);
	return 0;
}

