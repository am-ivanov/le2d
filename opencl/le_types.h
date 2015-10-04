#ifndef LE_TYPES_H
#define LE_TYPES_H

#ifdef USE_DOUBLE
	#pragma OPENCL EXTENSION cl_khr_fp64 : enable
	typedef double real;
#else
	typedef float real;
#endif

typedef int int_t;

#define NODE_SIZE 5
#define ind(i, j) ((i) + (j) * (nx))
#ifdef USE_AOS
	#define ind_vx(i, j) (*(grid + ind(i, j) * NODE_SIZE + 0))
	#define ind_vy(i, j) (*(grid + ind(i, j) * NODE_SIZE + 1))
	#define ind_sxx(i, j) (*(grid + ind(i, j) * NODE_SIZE + 2))
	#define ind_sxy(i, j) (*(grid + ind(i, j) * NODE_SIZE + 3))
	#define ind_syy(i, j) (*(grid + ind(i, j) * NODE_SIZE + 4))
#else // SOA
	#define ind_vx(i, j) (*(grid + ind(i, j) + 0 * nx * ny))
	#define ind_vy(i, j) (*(grid + ind(i, j)+ 1 * nx * ny))
	#define ind_sxx(i, j) (*(grid + ind(i, j) + 2 * nx * ny))
	#define ind_sxy(i, j) (*(grid + ind(i, j) + 3 * nx * ny))
	#define ind_syy(i, j) (*(grid + ind(i, j) + 4 * nx * ny))
#endif

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

#endif // LE_TYPES_H