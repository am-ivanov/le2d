/*
 * Author: Nikolay Khokhlov <k_h@inbox.ru>, (C) 2012
 * Compile: gcc -O2 -o le2d *.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include "gen.h"

//#include <omp.h>
#include "le_core.h"


#ifdef SAVE_EVERY_STEPS
const int save_every_step = 1;
#else
const int save_every_step = 0;
#endif

int main(int argc, char *argv[])
{
	if (argc != 4) {
		printf("Usage: %s nx ny steps.\n", argv[0]);
		return 1;
	}
	int i;
	le_point2 n = {atoi(argv[1]), atoi(argv[2])};
	int steps = atoi(argv[3]);
	le_task task;
	le_material mat;
	le_vec2 h = {1.0, 1.0};
	real dt = 0.3;
	real xxx = n.x / 2;
	real yyy = n.y / 2;
	le_vec2 center = {xxx, yyy};
	char name[1000];

	double t;
	unsigned long cc;

	//omp_set_num_threads(atoi(argv[4]));

	/*
	 * Init material.
	 */
	le_init_material(2.0, 1.0, 1.5, &mat);

	/*
	 * Init task.
	 */
	le_init_task(&task, dt, h, mat, n, ST_AOS);

	/*
	 * Initial condition.
	 */
	le_set_ball(&task, center, 10.0, 1.0);

	//int* indX = (int*) malloc(2 * task.n.x * task.n.y * steps * sizeof(int));
	//int* indY = (int*) malloc(2 * task.n.x * task.n.y * steps * sizeof(int));
	//int* indT = (int*) malloc(2 * task.n.x * task.n.y * steps * sizeof(int));

	//while (1) {

  // optimized
	memset(task.grid, 0, sizeof(le_node) * n.x * n.y);
	memset(task.grid2, 0, sizeof(le_node) * n.x * n.y);
	le_init_task(&task, dt, h, mat, n, ST_AOS);
	le_set_ball(&task, center, 10.0, 1.0);

	CachedGen generator(steps * 2, n.x / 2, n.y / 2);

	t = timer();
	le_step_rec(&task, steps, generator);
	t = timer() - t;
	printf("%d %d %d %f\n", n.x, n.y, steps, t);

	// original
	memset(task.grid, 0, sizeof(le_node) * n.x * n.y);
	memset(task.grid2, 0, sizeof(le_node) * n.x * n.y);
	le_init_task(&task, dt, h, mat, n, ST_AOS);
	le_set_ball(&task, center, 10.0, 1.0);

	t = timer();
	le_step_all2(&task, steps);
	t = timer() - t;
	printf("%d %d %d %f\n", n.x, n.y, steps, t);
	//memset(task.grid, 0, sizeof(le_node) * n.x * n.y);
	//memset(task.grid2, 0, sizeof(le_node) * n.x * n.y);
	//le_init_task(&task, dt, h, mat, n, ST_AOS);
	//le_set_ball(&task, center, 10.0, 1.0);

	//t = timer();
	//le_step_all(&task, steps);
	//t = timer() - t;
	//printf("%d %d %d %f\n", n.x, n.y, steps, t);

	//int k = 0;
	//for (int t = 0; t < 2 * steps; ++t) {
	//	for (int j = 0; j < task.n.y; ++j) {
	//		for (int i = 0; i < task.n.x; ++i) {
	//			indX[k] = i;
	//			indY[k] = j;
	//			indT[k] = t;
	//			++k;
	//		}
	//	}
	//}

	//memset(task.grid, 0, sizeof(le_node) * n.x * n.y);
	//memset(task.grid2, 0, sizeof(le_node) * n.x * n.y);
	//le_init_task(&task, dt, h, mat, n, ST_AOS);
	//le_set_ball(&task, center, 10.0, 1.0);

	//t = timer();
	//le_step_all3(&task, steps, indX, indY, indT);
	//t = timer() - t;
	//printf("%d %d %d %f\n", n.x, n.y, steps, t);

	//{
	//	le_task* t = &task;
	//	int Lx = 0;
	//	int Ux = t->n.x / 2;
	//	int Ly = 0;
	//	int Uy = t->n.y / 2;
	//	int Lt = 0;
	//	int Ut = steps * 2;
	//	int h = 16;

	//	int k = 0;
	//	for (int P = roundUp(3 * (Lt - h + 1), h) / h; P <= roundDown(3 * Ut, h) / h; P += 1) {
	//		int a_top = roundDown(3 * Ux + P*h + 2 * h - 2, 3 * h) / (3 * h);
	//		int a_bot = roundUp(3 * Lx + P*h - 2 * h + 2, 3 * h) / (3 * h);
	//		int b_top = roundDown(3 * Uy + P*h + 2 * h - 2, 3 * h) / (3 * h);
	//		int b_bot = roundUp(3 * Ly + P*h - 2 * h + 2, 3 * h) / (3 * h);
	//		for (int b = b_bot; b <= b_top; b += 1) {
	//			for (int a = a_bot; a <= a_top; a += 1) {
	//				int c = P - a - b;
	//				int tt_bot = MAX(roundUp((a + b + c) * h, 3), 3 * Lt) / 3;
	//				int tt_top = MIN(roundDown((a + b + c + 3) * h - 3, 3), 3 * (Ut - 1)) / 3;
	//				for (int tt = tt_bot; tt <= tt_top; tt += 1) {
	//					int j_bot = MAX(-tt + b * h, Ly);
	//					int j_top = MIN(-tt + (b + 1) * h - 1, Uy - 1);
	//					for (int j = j_bot; j <= j_top; ++j) {
	//						int i_bot = MAX(MAX(-tt + a * h, tt - j - (c + 1) * h + 1), Lx);
	//						int i_top = MIN(MIN(-tt + (a + 1) * h - 1, tt - j  - c * h), Ux - 1);
	//						for (int i = i_bot; i <= i_top; ++i) {
	//							for (int j1 = j * 2; j1 < j * 2 + 2; ++j1) {
	//								for (int i1 = i * 2; i1 < i * 2 + 2; ++i1) {
	//									indX[k] = i1;
	//									indY[k] = j1;
	//									indT[k] = tt;
	//									++k;
	//								}
	//							}
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}
	//}

	//memset(task.grid, 0, sizeof(le_node) * n.x * n.y);
	//memset(task.grid2, 0, sizeof(le_node) * n.x * n.y);
	//le_init_task(&task, dt, h, mat, n, ST_AOS);
	//le_set_ball(&task, center, 10.0, 1.0);

	//t = timer();
	//le_step_all3(&task, steps, indX, indY, indT);
	//t = timer() - t;
	//printf("%d %d %d %f\n", n.x, n.y, steps, t);
	//}

	//free(indX);
	//free(indY);
	//free(indT);

	/*
	 * Save last step.
	 */
	le_save_task(&task, "result.vtk");

	/*
	 * Free memory.
	 */
	le_free_task(&task);
	return 0;
}
