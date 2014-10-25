#ifndef LE_CORE_H
#define LE_CORE_H

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


/*
 * Storage types:
 * ST_SOA - Structure of srrays.
 * ST_AOS - Array of structures.
 */
#define ST_SOA 0
#define ST_AOS 1

/* Real type. */
#ifdef USE_DOUBLE
typedef double real;

#ifdef USE_SSE
	#define REAL_PER_SSE 2
#elif USE_AVX
	#define REAL_PER_SSE 4
#else
	#define REAL_PER_SSE 1
#endif

#else
typedef float real;

#ifdef USE_SSE
	#define REAL_PER_SSE 4
#elif USE_AVX
	#define REAL_PER_SSE 8
#else
	#define REAL_PER_SSE 1
#endif


#endif



/* No sse */
#if REAL_PER_SSE == 1
	#define sse_t real
	#define sse_add(x, y) (x + y)
	#define sse_sub(x, y) (x - (y))
	#define sse_mul(x, y) ((x) * (y))
	#define sse_div(x, y) ((x) / (y))
	#define sse_set(x) (x)
	#define sse_min(x, y) ((x) > (y) ? (y) : (x))
	#define sse_max(x, y) ((x) < (y) ? (y) : (x))
	#define sse_malloc(x) (malloc((x)))
#endif

/* double + sse */
#if REAL_PER_SSE == 2
	#define sse_t __m128d
	#define sse_add _mm_add_pd
	#define sse_sub _mm_sub_pd
	#define sse_mul _mm_mul_pd
	#define sse_div _mm_div_pd
	#define sse_set _mm_set1_pd
	#define sse_min _mm_min_pd
	#define sse_max _mm_max_pd
	#define sse_malloc(x) (memalign(16, (x)))
#endif

#if REAL_PER_SSE == 4
	/* double + avx */
	#ifdef USE_DOUBLE
		#define sse_t __m256d
		#define sse_add _mm256_add_pd
		#define sse_sub _mm256_sub_pd
		#define sse_mul _mm256_mul_pd
		#define sse_div _mm256_div_pd
		#define sse_set _mm256_set1_pd
		#define sse_min _mm256_min_pd
		#define sse_max _mm256_max_pd
		#define sse_malloc(x) (memalign(32, (x)))
	#else /* float + sse */
		#define sse_t __m128
		#define sse_add _mm_add_ps
		#define sse_sub _mm_sub_ps
		#define sse_mul _mm_mul_ps
		#define sse_div _mm_div_ps
		#define sse_set _mm_set1_ps
		#define sse_min _mm_min_ps
		#define sse_max _mm_max_ps
		#define sse_malloc(x) (memalign(16, (x)))
	#endif
#endif

/* float + avx */
#if REAL_PER_SSE == 8
	#define sse_t __m256
	#define sse_add _mm256_add_ps
	#define sse_sub _mm256_sub_ps
	#define sse_mul _mm256_mul_ps
	#define sse_div _mm256_div_ps
	#define sse_set _mm256_set1_ps
	#define sse_min _mm256_min_ps
	#define sse_max _mm256_max_ps
	#define sse_malloc(x) (memalign(32, (x)))
#endif


/* Integer type. */
typedef int int_t;

/* 2d vector struct. */
typedef struct {
	real x, y;
} le_vec2;

/* 
 * Symmetrical tensor 2d.
 * sxx sxy
 * sxy sxx
 */
typedef struct {
	real xx, xy, yy;
} le_smatrix2;

/* Integer 2d point. */
typedef struct {
	int_t x, y;
} le_point2;

/*
 * Elastic node structure.
 * Store velocity and stress tensor (http://en.wikipedia.org/wiki/Stress_%28mechanics%29).
 */
typedef struct {
	le_vec2 v;
	le_smatrix2 s;
} le_node;

/*
 * Riemann invariant structure (http://en.wikipedia.org/wiki/Riemann_invariant).
 */
typedef struct {
	real w1, w2, w3, w4/*, w5*/;
} le_w;

/*
 * Elastic material structure.
 * We use simple approach when all parameters are constant on whole region.
 * c1 - speed of P-wave http://en.wikipedia.org/wiki/P-wave.
 * c2 - speed of S-wave http://en.wikipedia.org/wiki/S-wave.
 * rho - density.
 */
typedef struct {
	real c1, c2, rho;
	
	/*
	 * Some cached values to speedup calculations.
	 */
	real irhoc1; // 1.0 / (c1 * rho)
	real irhoc2; // 1.0 / (c2 * rho)
	real rhoc1; // c1 * rho
	real rhoc2; // c2 * rho
	real rhoc3; // c3 * rho
} le_material;

/*
 * Structure for storing all parameters of task.
 */
typedef struct {
	/* Time step.*/
	real dt;
	
	/* Grid spacing. */
	le_vec2 h;
	
	/* Number of nodes ing grid on each axis. */
	le_point2 n;
	
	/* Material. */
	le_material mat;
	
	/*
	 * Storage type.
	 * Array of structure or structure of arrays.
	 */
	int stype;
	
	/* Grid data (nodes). */
	real *grid;
} le_task;

/*
 * Create material and init all fields of structure.
 */
void le_init_material(const real c1, const real c2, const real rho, le_material *m);

/* Create task with given parameters. Allocate memory for nodes. */
void le_init_task(le_task *task, const real dt, const le_vec2 h, const le_material mat, const le_point2 n, const int stype);

/* Free memory. */
void le_free_task(le_task *task);

/*
 * Set initial disturbance on the grid.
 */
void le_set_ball(le_task *t, const le_vec2 c, const real r, const real s);

/*
 * Save grid to legasy VTK format (http://www.vtk.org/VTK/img/file-formats.pdf).
 * You can use ParaView (http://www.paraview.org/),
 * MayaVi (http://mayavi.sourceforge.net/) or VisIt (https://wci.llnl.gov/codes/visit/)
 * to visualize results.
 * Return: 0 - all ok, 1 - error.
 */
int le_save_task(le_task *task, const char *file);

/*
 * One time step of difference scheme.
 */
void le_step(le_task *task);

/*
 * One time step of difference scheme.
 * Cache friendly version.
 */
void le_step_cf(le_task *task, const int cfs);

/*
 * One time step of difference scheme.
 * SOA version.
 */
void le_step_soa(le_task *task);

#endif //LE_CORE_H
