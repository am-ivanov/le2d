#ifndef LE_CORE_H
#define LE_CORE_H

#include "le_types.h"
#include <sys/time.h>

/*
 * Author: Nikolay Khokhlov <k_h@inbox.ru>, (C) 2012
 * Base include file.
 * All data structures and functions.
 * Prefix le_* used (Linear Elastic) for all structures and function.
 */

static __inline__ unsigned long getCC(void)
{
	unsigned a, d;
	asm volatile("rdtsc" : "=a" (a), "=d" (d));
	return ((unsigned long)a) | (((unsigned long)d) << 32);
}


static __inline__ double timer()
{
	struct timeval theStartTime;
	gettimeofday(&theStartTime, NULL);
	return theStartTime.tv_sec + 1e-6 * theStartTime.tv_usec;
}

/* Create material and init all fields of structure. */
void le_init_material(const real c1, const real c2, const real rho, le_material *m);

/* Create task with given parameters. Allocate memory for nodes. */
void le_init_task(le_task *task, const real dt, const le_vec2 h, const le_material mat, const le_point2 n);

/* Free memory. */
void le_free_task(le_task *task);

/* Set initial disturbance on the grid. */
void le_set_ball(le_task *t, const le_vec2 c, const real r, const real s);

/*
 * Save grid to legasy VTK format (http://www.vtk.org/VTK/img/file-formats.pdf).
 * You can use ParaView (http://www.paraview.org/),
 * MayaVi (http://mayavi.sourceforge.net/) or VisIt (https://wci.llnl.gov/codes/visit/)
 * to visualize results.
 * Return: 0 - all ok, 1 - error.
 */
int le_save_task(le_task *task, const char *file);

/* One time step of difference scheme. */
double le_step(le_task *task, int steps, int devices);

#endif //LE_CORE_H

