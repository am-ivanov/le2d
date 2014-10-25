#include <malloc.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <omp.h>

#include "le_core.h"

#define vnorm(v) (sqrt(v.x * v.x + v.y * v.y))

#define TVD2_EPS 1e-6

inline real le_min(real a, real b) { return a > b ? b : a; }
inline real le_max(real a, real b) { return a > b ? a : b; }
inline real le_max3(real a, real b, real c) { return le_max(a, le_max(b, c)); }

#define limiter_minmod(r) (le_max(0.0, le_min(1.0, (r))))
#define limiter_cir(r) (0.0)
#define limiter_superbee(r) (le_max3(0.0, le_min(1.0, 2.0 * r), le_min(2.0, r)))

#define limiter limiter_superbee

inline real tvd2(const real c, const real u_2, const real u_1, const real u, const real u1)
{
	const real eps = TVD2_EPS;
	real r1 = u - u_1;
	r1 += eps;
	real r2 = u1 - u;
	r2 += eps;
	const real r = r1 / r2;

	real r_1 = u_1 - u_2;
	r_1 += eps;
	r_1 /= r1;

	real f12 = r2 * limiter(r);
	real f_12 = r1 * limiter(r_1);

	const real k = 0.5 * (1 - c);

	return c * ((f_12 - f12) * k - r1);
}

void le_set_ball(le_task *t, const le_vec2 c, const real r, const real s)
{
	int i, j;
	for (i = 0; i < t->n.x; i++) {
		for (j = 0; j < t->n.y; j++) {
			le_vec2 x = {t->h.x * i, t->h.y * j};
			le_vec2 d = {x.x - c.x, x.y - c.y};
			if (vnorm(d) < r) {
				/* Set pressure disturbance */
				int_t nx = t->n.x;
				int_t ny = t->n.y;
				real* grid = t->grid;
				ind_sxx(i, j) = s;
				ind_syy(i, j) = s;
			}
		}
	}
}

/*
 * Write float to file and reverse byte order.
 */
void write_float(FILE* f, const float v)
{
	union {
		float f;
		unsigned char b[4];
	} dat1, dat2;
	dat1.f = v;
	dat2.b[0] = dat1.b[3];
	dat2.b[1] = dat1.b[2];
	dat2.b[2] = dat1.b[1];
	dat2.b[3] = dat1.b[0];
	fwrite(dat2.b, sizeof(unsigned char), 4, f);
}

void le_init_task(le_task *task, const real dt, const le_vec2 h, const le_material mat, const le_point2 n)
{
	task->dt  = dt;
	task->h   = h;
	task->mat = mat;
	task->n   = n;
	task->grid = (real*) malloc(sizeof(real) * NODE_SIZE * n.x * n.y);
	memset(task->grid, 0, NODE_SIZE * n.x * n.y);
}

void le_free_task(le_task* task)
{
	free(task->grid);
}

int le_save_task(le_task *t, const char *file)
{
	int i, j;
	FILE *fp = fopen(file, "w");
	if (fp == NULL) {
		perror("Failed to open file");
		return 1;
	}
	fprintf(fp, "# vtk DataFile Version 3.0\n");
	fprintf(fp, "Created by le_save_task\n");
	fprintf(fp, "BINARY\n");
	fprintf(fp, "DATASET STRUCTURED_POINTS\n");
	fprintf(fp, "DIMENSIONS %d %d 1\n", t->n.x, t->n.y);
	fprintf(fp, "SPACING %f %f 0.0\n", t->h.x, t->h.y);
	fprintf(fp, "ORIGIN 0.0 0.0 0.0\n");
	fprintf(fp, "POINT_DATA %d\n", t->n.x * t->n.y);

	/* velocity */
	fprintf(fp, "SCALARS v float 1\n");
	fprintf(fp, "LOOKUP_TABLE v_table\n");
	for (j = 0; j < t->n.y; j++) {
		for (i = 0; i < t->n.x; i++) {
			float v;
			const int_t nx = t->n.x;
			const int_t ny = t->n.y;
			const real* grid = t->grid;
			le_vec2 vt = { ind_vx(i, j), ind_vy(i, j) };
			v = vnorm(vt);
			write_float(fp, v);
		}
	}
	/*
	 * You can use the same code for saving other variables.
	 */
	fclose(fp);
	return 0;
}

void le_init_material(const real c1, const real c2, const real rho, le_material *m)
{
	m->c1 = c1;
	m->c2 = c2;
	m->rho = rho;

	/* Cached values. */
	m->irhoc1 = 1.0 / (c1 * rho);
	m->irhoc2 = 1.0 / (c2 * rho);
	m->rhoc1 = c1 * rho;
	m->rhoc2 = c2 * rho;
	real mu = rho * c2 * c2;
	real la = rho * c1 * c1 - 2.0 * mu;
	m->rhoc3 = rho * c1 * la / (la + 2.0 * mu);
}

inline void omega_x(const le_material *m, const real vx, const real vy, const real sxx, const real sxy, le_w *w)
{
	const real N00T = sxx * m->irhoc1;
	const real N01T = sxy * m->irhoc2;

	w->w1 = vx - N00T;
	w->w2 = vx + N00T;
	w->w3 = vy - N01T;
	w->w4 = vy + N01T;
}

inline void omega_y(const le_material *m, const real vx, const real vy, const real sxy, const real syy, le_w *w)
{
	const real N00T = syy * m->irhoc1;
	const real N01T = sxy * m->irhoc2;

	w->w1 = vy - N00T;
	w->w2 = vy + N00T;
	w->w3 = vx - N01T;
	w->w4 = vx + N01T;
}

inline void inc_x(const le_material *m, real* vx, real* vy, real* sxx, real* sxy, real* syy, const le_w *d)
{
	const real d1 = 0.5 * d->w1;
	const real d2 = 0.5 * d->w2;
	const real d3 = 0.5 * d->w3;
	const real d4 = 0.5 * d->w4;

	*vx += d1 + d2;
	*vy += d3 + d4;

	*syy += (d2 - d1) * m->rhoc3;
	*sxx += (d2 - d1) * m->rhoc1;
	*sxy += (d4 - d3) * m->rhoc2;
}

inline void inc_y(const le_material *m, real* vx, real* vy, real* sxx, real* sxy, real* syy, const le_w *d)
{
	const real d1 = 0.5 * d->w1;
	const real d2 = 0.5 * d->w2;
	const real d3 = 0.5 * d->w3;
	const real d4 = 0.5 * d->w4;

	*vy += d1 + d2;
	*vx += d3 + d4;

	*syy += (d2 - d1) * m->rhoc1;
	*sxx += (d2 - d1) * m->rhoc3;
	*sxy += (d4 - d3) * m->rhoc2;
}

inline void reconstruct(const le_w ppu, const le_w pu, const le_w u, const le_w nu, const le_w nnu, const real k1, const real k2, le_w *d)
{
	d->w1 = tvd2(k1, ppu.w1, pu.w1, u.w1, nu.w1); // c1
	d->w2 = tvd2(k1, nnu.w2, nu.w2, u.w2, pu.w2); // -c1
	d->w3 = tvd2(k2, ppu.w3, pu.w3, u.w3, nu.w3); // c2
	d->w4 = tvd2(k2, nnu.w4, nu.w4, u.w4, pu.w4); // -c2
}

void le_step_x(le_task *t)
{
	#pragma omp parallel
	{
	int i, j;

	const real k1 = t->dt * t->mat.c1 / t->h.x;
	const real k2 = t->dt * t->mat.c2 / t->h.x;
	le_material mat = t->mat;
	real* grid = t->grid;
	int_t nx = t->n.x;
	int_t ny = t->n.y;
	le_w d;

	real* w = (real*) malloc(sizeof(real) * W_SIZE * (nx+4));

	#define ind_w(i, k) (w[(k) * (nx + 4) + ((i)+2)])

	int parallel_width = nx / omp_get_num_threads();

	#pragma omp for schedule(static, parallel_width)
	for (j = 0; j < ny; j++) {
		#ifdef AUTOVECT
		#pragma omp simd
		#endif
		for (i = 0; i < nx; i++)
		{
			const real N00T = ind_sxx(i, j) * mat.irhoc1;
			const real N01T = ind_sxy(i, j) * mat.irhoc2;
			ind_w(i, 0) = ind_vx(i, j) - N00T;
			ind_w(i, 1) = ind_vx(i, j) + N00T;
			ind_w(i, 2) = ind_vy(i, j) - N01T;
			ind_w(i, 3) = ind_vy(i, j) + N01T;
		}
		ind_w(nx+1, 0) = ind_w(nx, 0) = ind_w(nx-1, 0);
		ind_w(nx+1, 1) = ind_w(nx, 1) = ind_w(nx-1, 1);
		ind_w(nx+1, 2) = ind_w(nx, 2) = ind_w(nx-1, 2);
		ind_w(nx+1, 3) = ind_w(nx, 3) = ind_w(nx-1, 3);
		ind_w(-2, 0) = ind_w(-1, 0) = ind_w(0, 0);
		ind_w(-2, 1) = ind_w(-1, 1) = ind_w(0, 1);
		ind_w(-2, 2) = ind_w(-1, 2) = ind_w(0, 2);
		ind_w(-2, 3) = ind_w(-1, 3) = ind_w(0, 3);
		#ifdef AUTOVECT
		#pragma omp simd
		#endif
		for (i = 0; i < nx; i++) {
			d.w1 = tvd2(k1, ind_w(i-2, 0), ind_w(i-1, 0), ind_w(i, 0), ind_w(i+1, 0));
			d.w2 = tvd2(k1, ind_w(i+2, 1), ind_w(i+1, 1), ind_w(i, 1), ind_w(i-1, 1));
			d.w3 = tvd2(k2, ind_w(i-2, 2), ind_w(i-1, 2), ind_w(i, 2), ind_w(i+1, 2));
			d.w4 = tvd2(k2, ind_w(i+2, 3), ind_w(i+1, 3), ind_w(i, 3), ind_w(i-1, 3));
			inc_x(&mat, ind_all(i, j), &d);
		}
	}
	free(w);
	}//end of omp parallel
#undef ind_w
}

void le_step_y(le_task *t)
{
	#define ind_w_2(i, k) (w_2[(k) * parallel_width + (i)])
	#define ind_w_1(i, k) (w_1[(k) * parallel_width + (i)])
	#define ind_w(i, k) (w[(k) * parallel_width + (i)])
	#define ind_w1(i, k) (w1[(k) * parallel_width + (i)])
	#define ind_w2(i, k) (w2[(k) * parallel_width + (i)])

	#define om_y(i, j, wn) \
	{\
		const real N00T = ind_syy(i + parallel_start, j) * mat.irhoc1;\
		const real N01T = ind_sxy(i + parallel_start, j) * mat.irhoc2;\
		ind_##wn(i, 0) = ind_vy(i + parallel_start, j) - N00T;\
		ind_##wn(i, 1) = ind_vy(i + parallel_start, j) + N00T;\
		ind_##wn(i, 2) = ind_vx(i + parallel_start, j) - N01T;\
		ind_##wn(i, 3) = ind_vx(i + parallel_start, j) + N01T;\
	}

	#define w_copy(from, to) \
	{\
		ind_##to(i, 0) = ind_##from(i, 0);\
		ind_##to(i, 1) = ind_##from(i, 1);\
		ind_##to(i, 2) = ind_##from(i, 2);\
		ind_##to(i, 3) = ind_##from(i, 3);\
	}

	#pragma omp parallel
	{
	int i, j;

	const real k1 = t->dt * t->mat.c1 / t->h.y;
	const real k2 = t->dt * t->mat.c2 / t->h.y;

	le_material mat = t->mat;
	real* grid = t->grid;
	int_t nx = t->n.x;
	int_t ny = t->n.y;
	le_w d;
	real* temp;

	int thread_num = omp_get_thread_num();
	int num_threads = omp_get_num_threads();
	int parallel_width = nx / num_threads;
	int parallel_start = thread_num * parallel_width;
	int parallel_end = parallel_start + parallel_width;
	if (thread_num == num_threads-1) {
		parallel_width = nx - parallel_start;
		parallel_end = nx;
	}

	real* w_2 = (real*) malloc(sizeof(real) * W_SIZE * parallel_width);
	real* w_1 = (real*) malloc(sizeof(real) * W_SIZE * parallel_width);
	real* w = (real*) malloc(sizeof(real) * W_SIZE * parallel_width);
	real* w1 = (real*) malloc(sizeof(real) * W_SIZE * parallel_width);
	real* w2 = (real*) malloc(sizeof(real) * W_SIZE * parallel_width);

	#ifdef AUTOVECT
	#pragma omp simd
	#endif
	for (i = 0; i < parallel_width; i++)
	{
		om_y(i, 0, w);
		om_y(i, 1, w1);
		om_y(i, 2, w2);
		w_copy(w, w_2);
		w_copy(w, w_1);
	}

	for (j = 0; j < ny; j++) {
		#ifdef AUTOVECT
		#pragma omp simd
		#endif
		for (i = 0; i < parallel_width; i++) {
			d.w1 = tvd2(k1, ind_w_2(i, 0), ind_w_1(i, 0), ind_w(i, 0), ind_w1(i, 0));
			d.w2 = tvd2(k1, ind_w2(i, 1), ind_w1(i, 1), ind_w(i, 1), ind_w_1(i, 1));
			d.w3 = tvd2(k2, ind_w_2(i, 2), ind_w_1(i, 2), ind_w(i, 2), ind_w1(i, 2));
			d.w4 = tvd2(k2, ind_w2(i, 3), ind_w1(i, 3), ind_w(i, 3), ind_w_1(i, 3));
			inc_y(&mat, ind_all(i + parallel_start, j), &d);
		}
		temp = w_2;
		w_2 = w_1;
		w_1 = w;
		w = w1;
		w1 = w2;
		w2 = temp;
		if (j < ny - 3) {
			#ifdef AUTOVECT
			#pragma omp simd
			#endif
			for (i = 0; i < parallel_width; i++) {
				om_y(i, j + 3, w2);
			}
		}
		else {
			#ifdef AUTOVECT
			#pragma omp simd
			#endif
			for (i = 0; i < parallel_width; i++) {
				w_copy(w1, w2);
			}
		}
	}
	free(w_2);
	free(w_1);
	free(w);
	free(w1);
	free(w2);
	}// end of omp parallel

	#undef ind_w_2
	#undef ind_w_1
	#undef ind_w
	#undef ind_w1
	#undef ind_w2
	#undef om_y
	#undef w_copy
}

void le_step_y_as_x(le_task *t)
{
	#define ind_w_2(i, k) (w_2[(k) * nx + (i)])
	#define ind_w_1(i, k) (w_1[(k) * nx + (i)])
	#define ind_w(i, k) (w[(k) * nx + (i)])
	#define ind_w1(i, k) (w1[(k) * nx + (i)])
	#define ind_w2(i, k) (w2[(k) * nx + (i)])

	#define ind_wl0(i, k) (ln4[4 * (thread_num) + 0][(k) * nx + (i)])
	#define ind_wl1(i, k) (ln4[4 * (thread_num) + 1][(k) * nx + (i)])
	#define ind_wl2(i, k) (ln4[4 * (thread_num) + 2][(k) * nx + (i)])
	#define ind_wl3(i, k) (ln4[4 * (thread_num) + 3][(k) * nx + (i)])

	#define ind_wl_last(i, k) (ln4[4 * (thread_num+1) + j - (parallel_width-5)][(k) * nx + (i)])

	#define om_y(i, j, wn) \
	{\
		const real N00T = ind_syy(i, j) * mat.irhoc1;\
		const real N01T = ind_sxy(i, j) * mat.irhoc2;\
		ind_##wn(i, 0) = ind_vy(i, j) - N00T;\
		ind_##wn(i, 1) = ind_vy(i, j) + N00T;\
		ind_##wn(i, 2) = ind_vx(i, j) - N01T;\
		ind_##wn(i, 3) = ind_vx(i, j) + N01T;\
	}

	#define w_copy(from, to) \
	{\
		ind_##to(i, 0) = ind_##from(i, 0);\
		ind_##to(i, 1) = ind_##from(i, 1);\
		ind_##to(i, 2) = ind_##from(i, 2);\
		ind_##to(i, 3) = ind_##from(i, 3);\
	}

	int max_threads = omp_get_max_threads();
	//printf("%d", max_threads);
	real** ln4 = malloc(sizeof(real*) * (max_threads+1) * 4); // 4 lines between threads

	#pragma omp parallel shared(ln4)
	{
	int i, j;

	const real k1 = t->dt * t->mat.c1 / t->h.y;
	const real k2 = t->dt * t->mat.c2 / t->h.y;

	le_material mat = t->mat;
	real* grid = t->grid;
	int_t nx = t->n.x;
	int_t ny = t->n.y;
	le_w d;
	real* temp;

	real* w_2 = (real*) malloc(sizeof(real) * W_SIZE * nx);
	real* w_1 = (real*) malloc(sizeof(real) * W_SIZE * nx);
	real* w = (real*) malloc(sizeof(real) * W_SIZE * nx);
	real* w1 = (real*) malloc(sizeof(real) * W_SIZE * nx);
	real* w2 = (real*) malloc(sizeof(real) * W_SIZE * nx);

	int thread_num = omp_get_thread_num();
	int num_threads = omp_get_num_threads();
	int parallel_width = ny / num_threads;
	int parallel_start = thread_num * parallel_width;
	int parallel_end = parallel_start + parallel_width;
	if (thread_num == num_threads-1) {
		parallel_width = ny - parallel_start;
		parallel_end = ny;
	}

	//printf("I am %d and i calc lines from %d to %d, width: %d\n", thread_num, parallel_start, parallel_end, parallel_width);

	ln4[4 * (thread_num+1) + 0] = (real*) malloc(sizeof(real) * W_SIZE * nx);
	ln4[4 * (thread_num+1) + 1] = (real*) malloc(sizeof(real) * W_SIZE * nx);
	ln4[4 * (thread_num+1) + 2] = (real*) malloc(sizeof(real) * W_SIZE * nx);
	ln4[4 * (thread_num+1) + 3] = (real*) malloc(sizeof(real) * W_SIZE * nx);

#pragma omp barrier

	if (thread_num != 0) {/* this 'if' is bad in some cases */
		#ifdef AUTOVECT
		#pragma omp simd
		#endif
		for (i = 0; i < nx; i++) {
			om_y(i, parallel_start-2, w_2);
			w_copy(w_2, wl0);
			om_y(i, parallel_start-1, w_1);
			w_copy(w_1, wl1);
			om_y(i, parallel_start+0, w);
			w_copy(w, wl2);
			om_y(i, parallel_start+1, w1);
			w_copy(w1, wl3);
			om_y(i, parallel_start+2, w2);
		}
		//printf("I am %d and i save lines from %d\n", thread_num, &(ln4[4 * thread_num]));
	} else {
		#ifdef AUTOVECT
		#pragma omp simd
		#endif
		for (i = 0; i < nx; i++) {
			om_y(i, 0, w);
			om_y(i, 1, w1);
			om_y(i, 2, w2);
			w_copy(w, w_2);
			w_copy(w, w_1);
		}
	}

	for (j = 0; j < parallel_width-5; j++) {
		#ifdef AUTOVECT
		#pragma omp simd
		#endif
		for (i = 0; i < nx; i++) {
			d.w1 = tvd2(k1, ind_w_2(i, 0), ind_w_1(i, 0), ind_w(i, 0), ind_w1(i, 0));
			d.w2 = tvd2(k1, ind_w2(i, 1), ind_w1(i, 1), ind_w(i, 1), ind_w_1(i, 1));
			d.w3 = tvd2(k2, ind_w_2(i, 2), ind_w_1(i, 2), ind_w(i, 2), ind_w1(i, 2));
			d.w4 = tvd2(k2, ind_w2(i, 3), ind_w1(i, 3), ind_w(i, 3), ind_w_1(i, 3));
			inc_y(&mat, ind_all(i, j+parallel_start), &d);
		}
		temp = w_2;
		w_2 = w_1;
		w_1 = w;
		w = w1;
		w1 = w2;
		w2 = temp;
		#ifdef AUTOVECT
		#pragma omp simd
		#endif
		for (i = 0; i < nx; i++) {
			om_y(i, j+parallel_start+3, w2);
		}
	}

	#pragma omp barrier

	for (j = parallel_width-5; j < parallel_width; j++) {
		#ifdef AUTOVECT
		#pragma omp simd
		#endif
		for (i = 0; i < nx; i++) {
			d.w1 = tvd2(k1, ind_w_2(i, 0), ind_w_1(i, 0), ind_w(i, 0), ind_w1(i, 0));
			d.w2 = tvd2(k1, ind_w2(i, 1), ind_w1(i, 1), ind_w(i, 1), ind_w_1(i, 1));
			d.w3 = tvd2(k2, ind_w_2(i, 2), ind_w_1(i, 2), ind_w(i, 2), ind_w1(i, 2));
			d.w4 = tvd2(k2, ind_w2(i, 3), ind_w1(i, 3), ind_w(i, 3), ind_w_1(i, 3));
			inc_y(&mat, ind_all(i, j+parallel_start), &d);
		}
		if (j != parallel_width - 1) {
			temp = w_2;
			w_2 = w_1;
			w_1 = w;
			w = w1;
			w1 = w2;
			w2 = temp;
			if (thread_num != num_threads - 1) { // if NOT last
				#ifdef AUTOVECT
				#pragma omp simd
				#endif
				for (i = 0; i < nx; i++) {
					w_copy(wl_last, w2);
					//w_copy(w2, wl_last);
				}
				//if (j == parallel_width-5) printf("I am %d and i use lines from %d\n", thread_num, &(ln4[4 * (thread_num+1)]));
			} else { // if last
				if (j < parallel_width - 3) {
					#ifdef AUTOVECT
					#pragma omp simd
					#endif
					for (i = 0; i < nx; i++) {
						om_y(i, j+parallel_start+3, w2);
					}
				} else {
					#ifdef AUTOVECT
					#pragma omp simd
					#endif
					for (i = 0; i < nx; i++) {
						w_copy(w1, w2);
					}
				}
			}
		}
	}

	free(w_2);
	free(w_1);
	free(w);
	free(w1);
	free(w2);

	free(ln4[4 * (thread_num+1) + 0]);
	free(ln4[4 * (thread_num+1) + 1]);
	free(ln4[4 * (thread_num+1) + 2]);
	free(ln4[4 * (thread_num+1) + 3]);

	} /* end of omp parallel */

	free(ln4);

	#undef ind_w_2
	#undef ind_w_1
	#undef ind_w
	#undef ind_w1
	#undef ind_w2
	#undef om_y
	#undef w_copy
}

void le_step_y_base(le_task *t)
{
#pragma omp parallel
	{
	int i, j;

	const real k1 = t->dt * t->mat.c1 / t->h.y;
	const real k2 = t->dt * t->mat.c2 / t->h.y;
	le_material mat = t->mat;
	real* grid = t->grid;
	int_t nx = t->n.x;
	int_t ny = t->n.y;
	le_w d;

	real* w = (real*) malloc(sizeof(real) * W_SIZE * (ny+4));

	#define ind_w(i, k) (w[(k) * (ny + 4) + ((i)+2)])

	#pragma omp for
	for (i = 0; i < nx; i++) {
		#ifdef AUTOVECT
		#pragma omp simd
		#endif
		for (j = 0; j < ny; j++) {
			const real N00T = ind_syy(i, j) * mat.irhoc1;
			const real N01T = ind_sxy(i, j) * mat.irhoc2;
			ind_w(j, 0) = ind_vy(i, j) - N00T;
			ind_w(j, 1) = ind_vy(i, j) + N00T;
			ind_w(j, 2) = ind_vx(i, j) - N01T;
			ind_w(j, 3) = ind_vx(i, j) + N01T;
		}
		ind_w(ny+1, 0) = ind_w(ny, 0) = ind_w(ny-1, 0);
		ind_w(ny+1, 1) = ind_w(ny, 1) = ind_w(ny-1, 1);
		ind_w(ny+1, 2) = ind_w(ny, 2) = ind_w(ny-1, 2);
		ind_w(ny+1, 3) = ind_w(ny, 3) = ind_w(ny-1, 3);
		ind_w(-2, 0) = ind_w(-1, 0) = ind_w(0, 0);
		ind_w(-2, 1) = ind_w(-1, 1) = ind_w(0, 1);
		ind_w(-2, 2) = ind_w(-1, 2) = ind_w(0, 2);
		ind_w(-2, 3) = ind_w(-1, 3) = ind_w(0, 3);
		#ifdef AUTOVECT
		#pragma omp simd
		#endif
		for (j = 0; j < ny; j++) {
			d.w1 = tvd2(k1, ind_w(j-2, 0), ind_w(j-1, 0), ind_w(j, 0), ind_w(j+1, 0));
			d.w2 = tvd2(k1, ind_w(j+2, 1), ind_w(j+1, 1), ind_w(j, 1), ind_w(j-1, 1));
			d.w3 = tvd2(k2, ind_w(j-2, 2), ind_w(j-1, 2), ind_w(j, 2), ind_w(j+1, 2));
			d.w4 = tvd2(k2, ind_w(j+2, 3), ind_w(j+1, 3), ind_w(j, 3), ind_w(j-1, 3));
			inc_y(&mat, ind_all(i, j), &d);
		}
	}
#undef ind_w

	}
}

void le_step(le_task *task)
{
	le_step_x(task);

#ifdef BASE
	le_step_y_base(task);
#elif Y_AS_X //parallelism
	le_step_y_as_x(task);
#else // Y_AS_Y parallelism
	le_step_y(task);
#endif
}
