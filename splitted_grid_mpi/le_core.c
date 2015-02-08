#include <malloc.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <mpi.h>

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
}

void le_init_grid(le_task *task)
{
	const le_point2 n = task->n;	
	task->grid = (real*) malloc(sizeof(real) * NODE_SIZE * n.x * n.y);
	memset(task->grid, 0, NODE_SIZE * n.x * n.y);
}

void le_free_task(le_task* task)
{
	free(task->grid);
}

void scatter_work(le_task* task, int size, int rank)
{
	// WARNING: MPI_DOUBLE works only with #define USE_DOUBLE
	int i;
	int div = task->n.y / size;
	int mod = task->n.y % size;
	task->pn.x = task->n.x;
	task->pn.y = div;
	if (rank < mod) 
		++task->pn.y;
	// allocate momry for working buffer
	task->pgrid = (real*) malloc(sizeof(real) * task->pn.x * (task->pn.y + 4) * NODE_SIZE); // additional 2 lines in beginning and 2 lines at the end
	// call MPI_Scatterv
	int* sendcounts = (int*) malloc(sizeof(int) * size);
	int* displs = (int*) malloc(sizeof(int) * size);
	displs[0] = 0;
	for (i = 0; i != mod; ++i)
		sendcounts[i] = task->pn.x * (div + 1) * NODE_SIZE;
	for (i = mod; i != size; ++i) 
		sendcounts[i] = task->pn.x * div * NODE_SIZE;
	displs[0] = 0;
	for (i = 1; i != size; ++i)
		displs[i] = sendcounts[i - 1] + displs[i - 1];
	MPI_Scatterv(task->grid, sendcounts, displs, MPI_DOUBLE,
		task->pgrid + task->pn.x * 2 * NODE_SIZE, task->pn.x * task->pn.y * NODE_SIZE, MPI_DOUBLE,
		0, MPI_COMM_WORLD);
	free(sendcounts);
	free(displs);
	// additional memory(5 additional lines)
	task->w = (real*) malloc(sizeof(real) * NODE_SIZE * task->n.x * 5);
}

void gather_work(le_task* task, int size)
{
	// WARNING: MPI_DOUBLE works only with #define USE_DOUBLE
	int i;
	int div = task->n.y / size;
	int mod = task->n.y % size;
	int* recvcounts = (int*) malloc(sizeof(int) * size);
	int* displs = (int*) malloc(sizeof(int) * size);
	displs[0] = 0;
	for (i = 0; i != mod; ++i)
		recvcounts[i] = task->pn.x * (div + 1) * NODE_SIZE;
	for (i = mod; i != size; ++i) 
		recvcounts[i] = task->pn.x * div * NODE_SIZE;
	displs[0] = 0;
	for (i = 1; i != size; ++i)
		displs[i] = recvcounts[i - 1] + displs[i - 1];
	
	MPI_Gatherv(task->pgrid + task->pn.x * 2 * NODE_SIZE, task->pn.x * task->pn.y * NODE_SIZE, MPI_DOUBLE,
		task->grid, recvcounts, displs, MPI_DOUBLE,
		0, MPI_COMM_WORLD);
	
	free(recvcounts);
	free(displs);	    
	
	free(task->pgrid);
	free(task->w);
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
	#define ind_w(i, k) (w[(k) * (nx + 4) + ((i)+2)])
	
	int i, j;
	
	const real k1 = t->dt * t->mat.c1 / t->h.x;
	const real k2 = t->dt * t->mat.c2 / t->h.x;
	le_material mat = t->mat;
	real* pgrid = t->pgrid;
	real* w = t->w;
	int_t nx = t->pn.x;
	int_t ny = t->pn.y;
	
	for (j = 2; j < ny + 2; j++) {
		#ifdef AUTOVECT
		#pragma omp simd
		#endif
		for (i = 0; i < nx; i++)
		{
			const real N00T = inds_sxx(i, j) * mat.irhoc1;
			const real N01T = inds_sxy(i, j) * mat.irhoc2;
			ind_w(i, 0) = inds_vx(i, j) - N00T;
			ind_w(i, 1) = inds_vx(i, j) + N00T;
			ind_w(i, 2) = inds_vy(i, j) - N01T;
			ind_w(i, 3) = inds_vy(i, j) + N01T;
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
			le_w d;
			d.w1 = tvd2(k1, ind_w(i-2, 0), ind_w(i-1, 0), ind_w(i, 0), ind_w(i+1, 0));
			d.w2 = tvd2(k1, ind_w(i+2, 1), ind_w(i+1, 1), ind_w(i, 1), ind_w(i-1, 1));
			d.w3 = tvd2(k2, ind_w(i-2, 2), ind_w(i-1, 2), ind_w(i, 2), ind_w(i+1, 2));
			d.w4 = tvd2(k2, ind_w(i+2, 3), ind_w(i+1, 3), ind_w(i, 3), ind_w(i-1, 3));
			inc_x(&mat, &inds_vx(i,j), &inds_vy(i,j), &inds_sxx(i,j), &inds_sxy(i,j), &inds_syy(i,j), &d);
		}
	}
	#undef ind_w
}

void le_step_y(le_task *t)
{
	#define ind_w(i, k, wn) ((wn)[(k) * nx + (i)])

	#define om_y(i, j, wn) {\
		const real N00T = inds_syy(i, j) * mat.irhoc1;\
		const real N01T = inds_sxy(i, j) * mat.irhoc2;\
		ind_w(i, 0, wn) = inds_vy(i, j) - N00T;\
		ind_w(i, 1, wn) = inds_vy(i, j) + N00T;\
		ind_w(i, 2, wn) = inds_vx(i, j) - N01T;\
		ind_w(i, 3, wn) = inds_vx(i, j) + N01T;\
	}
	
	#define om_border(i, wn, tconn) {\
		const real N00T = tconn[i * NODE_SIZE + 4] * mat.irhoc1;\
		const real N01T = tconn[i * NODE_SIZE + 3] * mat.irhoc2;\
		ind_w(i, 0, wn) = tconn[i * NODE_SIZE + 1] - N00T;\
		ind_w(i, 1, wn) = tconn[i * NODE_SIZE + 1] + N00T;\
		ind_w(i, 2, wn) = tconn[i * NODE_SIZE + 0] - N01T;\
		ind_w(i, 3, wn) = tconn[i * NODE_SIZE + 0] + N01T;\
	}

	#define w_copy(from, to) {\
		ind_w(i, 0, to) = ind_w(i, 0, from);\
		ind_w(i, 1, to) = ind_w(i, 1, from);\
		ind_w(i, 2, to) = ind_w(i, 2, from);\
		ind_w(i, 3, to) = ind_w(i, 3, from);\
	}
	
	#define split_save(j, tconn) {\
		tconn[i * NODE_SIZE + 0] = inds_vx(i, j);\
		tconn[i * NODE_SIZE + 1] = inds_vy(i, j);\
		tconn[i * NODE_SIZE + 2] = inds_sxx(i, j);\
		tconn[i * NODE_SIZE + 3] = inds_sxy(i, j);\
		tconn[i * NODE_SIZE + 4] = inds_syy(i, j);\
	}

	int i, j;

	const real k1 = t->dt * t->mat.c1 / t->h.y;
	const real k2 = t->dt * t->mat.c2 / t->h.y;

	le_material mat = t->mat;
	real* pgrid = t->pgrid;
	real* tw = t->w;
	int_t nx = t->pn.x;
	int_t ny = t->pn.y;

	real* w_2 = tw + 0 * NODE_SIZE * nx;
	real* w_1 = tw + 1 * NODE_SIZE * nx;
	real* w = tw + 2 * NODE_SIZE * nx;
	real* w1 = tw + 3 * NODE_SIZE * nx;
	real* w2 = tw + 4 * NODE_SIZE * nx;
	
	int size, rank, prev, next;
	MPI_Request request[4];
	MPI_Status status[4];
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	// send/recv first 2 lines and last 2 lines to previous and next process
	if (0 == rank) {
		// copy first additional 2 lines from first line
		memcpy(pgrid, pgrid + nx * NODE_SIZE * 2, nx * NODE_SIZE);
		memcpy(pgrid + nx * NODE_SIZE, pgrid + nx * NODE_SIZE * 2, nx * NODE_SIZE);
	}
	else {
		// send first 2 lines to prev process
		MPI_Isend(pgrid + nx * NODE_SIZE * 2, nx * NODE_SIZE * 2, MPI_DOUBLE, rank - 1, 5, MPI_COMM_WORLD, &request[0]);
		// recv first additional 2 lines from prev process
		MPI_Irecv(pgrid, nx * NODE_SIZE * 2, MPI_DOUBLE, rank - 1, 6, MPI_COMM_WORLD, &request[1]);
	}
	
	if (size - 1 == rank) {
		// copy last additional 2 lines from last line
		memcpy(pgrid + nx * NODE_SIZE * (2 + ny), pgrid + nx * NODE_SIZE * (1 + ny), nx * NODE_SIZE);
		memcpy(pgrid + nx * NODE_SIZE * (3 + ny), pgrid + nx * NODE_SIZE * (1 + ny), nx * NODE_SIZE);
	} 
	else {

		// send last 2 lines to next proess
		MPI_Isend(pgrid + nx * NODE_SIZE * (0 + ny), nx * NODE_SIZE * 2, MPI_DOUBLE, rank + 1, 6, MPI_COMM_WORLD, &request[2]);
		// recv last additional 2 lines from next process
		MPI_Irecv(pgrid + nx * NODE_SIZE * (2 + ny), nx * NODE_SIZE * 2, MPI_DOUBLE, rank + 1, 5, MPI_COMM_WORLD, &request[3]);
	}
	if (0 != rank)
		MPI_Waitall(2, request, status);
	if (size - 1 != rank)
		MPI_Waitall(2, request + 2, status + 2);

	for (i = 0; i < nx; i++)
	{
		om_y(i, 0, w_2);
		om_y(i, 1, w_1);
		om_y(i, 2, w);
		om_y(i, 3, w1);
	}

	for (j = 2; j < ny + 2; j++) {
		real* temp;
		for (i = 0; i < nx; i++) {
			om_y(i, j + 2, w2);
		}
		for (i = 0; i < nx; i++) {
			le_w d;
			d.w1 = tvd2(k1, ind_w(i, 0, w_2), ind_w(i, 0, w_1), ind_w(i, 0, w), ind_w(i, 0, w1));
			d.w2 = tvd2(k1, ind_w(i, 1, w2), ind_w(i, 1, w1), ind_w(i, 1, w), ind_w(i, 1, w_1));
			d.w3 = tvd2(k2, ind_w(i, 2, w_2), ind_w(i, 2, w_1), ind_w(i, 2, w), ind_w(i, 2, w1));
			d.w4 = tvd2(k2, ind_w(i, 3, w2), ind_w(i, 3, w1), ind_w(i, 3, w), ind_w(i, 3, w_1));
			inc_y(&mat, &inds_vx((i),(j)), &inds_vy((i),(j)), &inds_sxx((i),(j)), &inds_sxy((i),(j)), &inds_syy((i),(j)), &d);
		}
		temp = w_2;
		w_2 = w_1;
		w_1 = w;
		w = w1;
		w1 = w2;
		w2 = temp;
	}

	#undef ind_w
	#undef om_y
	#undef w_copy
}



void le_step(le_task *task)
{
	le_step_x(task);
	le_step_y(task);
}
