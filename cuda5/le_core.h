#ifndef LE_CORE_H
#define LE_CORE_H

#include <sys/time.h>

/*
 * Author: Nikolay Khokhlov <k_h@inbox.ru>, (C) 2012
 * Base include file.
 * All data structures and functions.
 * Prefix le_* used (Linear Elastic) for all structures and function.
 */

#ifdef DEBUG

#define checkCudaErrors(param) { cudaAssert((param), __FILE__, __LINE__); }
inline void cudaAssert(cudaError_t code, const char *file, int line)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"cudaAssert: %s %s %d\n", cudaGetErrorString(code), file, line);
      exit(code);
   }
}

#define printf_debug printf

#else

#define printf_debug
#define checkCudaErrors(param) { param; }

#endif

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

#ifdef USE_DOUBLE
	typedef double real;
#else
	typedef float real;
#endif

typedef int int_t;

#define NODE_SIZE 5
#define ind(i, j) ((i) + (j) * (nx))

#define ind_vx(i, j) (*(grid + ind(i, j) + 0 * nx * ny))
#define ind_vy(i, j) (*(grid + ind(i, j) + 1 * nx * ny))
#define ind_sxx(i, j) (*(grid + ind(i, j) + 2 * nx * ny))
#define ind_sxy(i, j) (*(grid + ind(i, j) + 3 * nx * ny))
#define ind_syy(i, j) (*(grid + ind(i, j) + 4 * nx * ny))

#define W_SIZE 4 /* w1, w2, w3, w4 */

typedef struct {
	real x, y;
} le_vec2;

typedef struct {
	int_t x, y;
} le_point2;

typedef struct {
	real w1, w2, w3, w4;
} le_w;

typedef struct {
	real c1, c2, rho;

	/* Some cached values to speedup calculations. */
	real irhoc1; // 1.0 / (c1 * rho)
	real irhoc2; // 1.0 / (c2 * rho)
	real rhoc1; // c1 * rho
	real rhoc2; // c2 * rho
	real rhoc3; // c3 * rho
} le_material;

/* Structure for storing all parameters of task. */
typedef struct {
	/* Time step.*/
	real dt;

	/* Grid spacing. */
	le_vec2 h;

	/* Number of nodes ing grid on each axis. */
	le_point2 n;

	/* Material. */
	le_material mat;

	/* Grid data (nodes). */
	real *grid;
} le_task;

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
double le_step(le_task *task, int steps);

#endif //LE_CORE_H


