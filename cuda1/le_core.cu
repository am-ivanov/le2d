#include <malloc.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include <cuda.h>
#include <cuda_runtime.h>

#include "le_core.h"

#ifdef DEBUG

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

#else

#define gpuErrchk(ans) { (ans); }

#endif

#define vnorm(v) (sqrt(v.x * v.x + v.y * v.y))

#define TVD2_EPS 1e-6

__device__ inline real le_min(real a, real b) { return a > b ? b : a; }
__device__ inline real le_max(real a, real b) { return a > b ? a : b; }
__device__ inline real le_max3(real a, real b, real c) { return le_max(a, le_max(b, c)); }

#define limiter_minmod(r) (le_max(0.0, le_min(1.0, (r))))
#define limiter_cir(r) (0.0)
#define limiter_superbee(r) (le_max3(0.0, le_min(1.0, 2.0 * r), le_min(2.0, r)))

#define limiter limiter_superbee

__device__ inline real tvd2(const real c, const real u_2, const real u_1, const real u, const real u1)
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
				//int_t ny = t->n.y;
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
	memset(task->grid, 0, sizeof(real) * NODE_SIZE * n.x * n.y);
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
			//const int_t ny = t->n.y;
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

__device__ inline void omega_x(const le_material *m, const real vx, const real vy, const real sxx, const real sxy, le_w *w)
{
	const real N00T = sxx * m->irhoc1;
	const real N01T = sxy * m->irhoc2;

	w->w1 = vx - N00T;
	w->w2 = vx + N00T;
	w->w3 = vy - N01T;
	w->w4 = vy + N01T;
}

__device__ inline void omega_y(const le_material *m, const real vx, const real vy, const real sxy, const real syy, le_w *w)
{
	const real N00T = syy * m->irhoc1;
	const real N01T = sxy * m->irhoc2;

	w->w1 = vy - N00T;
	w->w2 = vy + N00T;
	w->w3 = vx - N01T;
	w->w4 = vx + N01T;
}

__device__ inline void inc_x(const le_material *m, const real vx, const real vy, const real sxx, const real sxy, const real syy,
		real* new_vx, real* new_vy, real* new_sxx, real* new_sxy, real* new_syy, const le_w *d)
{
	const real d1 = 0.5 * d->w1;
	const real d2 = 0.5 * d->w2;
	const real d3 = 0.5 * d->w3;
	const real d4 = 0.5 * d->w4;

	*new_vx = vx + d1 + d2;
	*new_vy = vy + d3 + d4;

	*new_syy = syy + (d2 - d1) * m->rhoc3;
	*new_sxx = sxx + (d2 - d1) * m->rhoc1;
	*new_sxy = sxy + (d4 - d3) * m->rhoc2;
}

__device__ inline void inc_y(const le_material *m, const real vx, const real vy, const real sxx, const real sxy, const real syy,
		real* new_vx, real* new_vy, real* new_sxx, real* new_sxy, real* new_syy, const le_w *d)
{
	const real d1 = 0.5 * d->w1;
	const real d2 = 0.5 * d->w2;
	const real d3 = 0.5 * d->w3;
	const real d4 = 0.5 * d->w4;

	*new_vy = vy + d1 + d2;
	*new_vx = vx + d3 + d4;

	*new_syy = syy + (d2 - d1) * m->rhoc1;
	*new_sxx = sxx + (d2 - d1) * m->rhoc3;
	*new_sxy = sxy + (d4 - d3) * m->rhoc2;
}

__device__ inline void write_new_values(real* vx, real* vy, real* sxx, real* sxy, real* syy,
		const real new_vx, const real new_vy, const real new_sxx, const real new_sxy, const real new_syy)
{
	*vx = new_vx;
	*vy = new_vy;
	*sxx = new_sxx;
	*sxy = new_sxy;
	*syy = new_syy;
}

__device__ inline void reconstruct(const le_w ppu, const le_w pu, const le_w u, const le_w nu, const le_w nnu, const real k1, const real k2, le_w *d)
{
	d->w1 = tvd2(k1, ppu.w1, pu.w1, u.w1, nu.w1); // c1
	d->w2 = tvd2(k1, nnu.w2, nu.w2, u.w2, pu.w2); // -c1
	d->w3 = tvd2(k2, ppu.w3, pu.w3, u.w3, nu.w3); // c2
	d->w4 = tvd2(k2, nnu.w4, nu.w4, u.w4, pu.w4); // -c2
}

__global__ void le_step_x(le_task *t, real* in_grid, real* out_grid)
{
	int i, j;
	const real k1 = t->dt * t->mat.c1 / t->h.x;
	const real k2 = t->dt * t->mat.c2 / t->h.x;
	int nx = t->n.x;
	int ny = t->n.y;
	real* grid = in_grid;

	real vx, vy, sxx, sxy, syy;
	le_w w_2, w_1, w, w1, w2, d;

	i = threadIdx.x + blockDim.x * blockIdx.x;
	j = threadIdx.y + blockDim.y * blockIdx.y;

	if (i >= nx || j >= ny) return;

	if (2 <= i && i < nx - 2) {
		omega_x(&t->mat, ind_vx(i-2, j), ind_vy(i-2, j), ind_sxx(i-2, j), ind_sxy(i-2, j), &w_2);
		omega_x(&t->mat, ind_vx(i-1, j), ind_vy(i-1, j), ind_sxx(i-1, j), ind_sxy(i-1, j), &w_1);
		omega_x(&t->mat, ind_vx(i, j), ind_vy(i, j), ind_sxx(i, j), ind_sxy(i, j), &w);
		omega_x(&t->mat, ind_vx(i+1, j), ind_vy(i+1, j), ind_sxx(i+1, j), ind_sxy(i+1, j), &w1);
		omega_x(&t->mat, ind_vx(i+2, j), ind_vy(i+2, j), ind_sxx(i+2, j), ind_sxy(i+2, j), &w2);
	} else if (i == 0) {
		omega_x(&t->mat, ind_vx(i, j), ind_vy(i, j), ind_sxx(i, j), ind_sxy(i, j), &w);
		omega_x(&t->mat, ind_vx(i+1, j), ind_vy(i+1, j), ind_sxx(i+1, j), ind_sxy(i+1, j), &w1);
		omega_x(&t->mat, ind_vx(i+2, j), ind_vy(i+2, j), ind_sxx(i+2, j), ind_sxy(i+2, j), &w2);
		w_2 = w_1 = w;
	} else if (i == 1) {
		omega_x(&t->mat, ind_vx(i-1, j), ind_vy(i-1, j), ind_sxx(i-1, j), ind_sxy(i-1, j), &w_1);
		omega_x(&t->mat, ind_vx(i, j), ind_vy(i, j), ind_sxx(i, j), ind_sxy(i, j), &w);
		omega_x(&t->mat, ind_vx(i+1, j), ind_vy(i+1, j), ind_sxx(i+1, j), ind_sxy(i+1, j), &w1);
		omega_x(&t->mat, ind_vx(i+2, j), ind_vy(i+2, j), ind_sxx(i+2, j), ind_sxy(i+2, j), &w2);
		w_2 = w_1;
	} else if (i == nx - 2) {
		omega_x(&t->mat, ind_vx(i-2, j), ind_vy(i-2, j), ind_sxx(i-2, j), ind_sxy(i-2, j), &w_2);
		omega_x(&t->mat, ind_vx(i-1, j), ind_vy(i-1, j), ind_sxx(i-1, j), ind_sxy(i-1, j), &w_1);
		omega_x(&t->mat, ind_vx(i, j), ind_vy(i, j), ind_sxx(i, j), ind_sxy(i, j), &w);
		omega_x(&t->mat, ind_vx(i+1, j), ind_vy(i+1, j), ind_sxx(i+1, j), ind_sxy(i+1, j), &w1);
		w2 = w1;
	} else if (i == nx - 1) {
		omega_x(&t->mat, ind_vx(i-2, j), ind_vy(i-2, j), ind_sxx(i-2, j), ind_sxy(i-2, j), &w_2);
		omega_x(&t->mat, ind_vx(i-1, j), ind_vy(i-1, j), ind_sxx(i-1, j), ind_sxy(i-1, j), &w_1);
		omega_x(&t->mat, ind_vx(i, j), ind_vy(i, j), ind_sxx(i, j), ind_sxy(i, j), &w);
		w2 = w1 = w;
	}

	reconstruct(w_2, w_1, w, w1, w2, k1, k2, &d);
	inc_x(&t->mat, ind_vx(i,j), ind_vy(i,j), ind_sxx(i,j), ind_sxy(i,j), ind_syy(i,j), &vx, &vy, &sxx, &sxy, &syy, &d);
	grid = out_grid;
	write_new_values(&ind_vx(i,j), &ind_vy(i,j), &ind_sxx(i,j), &ind_sxy(i,j), &ind_syy(i,j), vx, vy, sxx, sxy, syy);
}

__global__ void le_step_y(le_task *t, real* in_grid, real* out_grid)
{
	int i, j;
	const real k1 = t->dt * t->mat.c1 / t->h.y;
	const real k2 = t->dt * t->mat.c2 / t->h.y;
	int nx = t->n.x;
	int ny = t->n.y;
	real* grid = in_grid;

	real vx, vy, sxx, sxy, syy;
	le_w w_2, w_1, w, w1, w2, d;

	i = threadIdx.x + blockDim.x * blockIdx.x;
	j = threadIdx.y + blockDim.y * blockIdx.y;

	if (i >= nx || j >= ny) return;

	if (2 <= j && j < ny - 2) {
		omega_y(&t->mat, ind_vx(i, j-2), ind_vy(i, j-2), ind_sxy(i, j-2), ind_syy(i, j-2), &w_2);
		omega_y(&t->mat, ind_vx(i, j-1), ind_vy(i, j-1), ind_sxy(i, j-1), ind_syy(i, j-1), &w_1);
		omega_y(&t->mat, ind_vx(i, j), ind_vy(i, j), ind_sxy(i, j), ind_syy(i, j), &w);
		omega_y(&t->mat, ind_vx(i, j+1), ind_vy(i, j+1), ind_sxy(i, j+1), ind_syy(i, j+1), &w1);
		omega_y(&t->mat, ind_vx(i, j+2), ind_vy(i, j+2), ind_sxy(i, j+2), ind_syy(i, j+2), &w2);
	} else if (j == 0) {
		omega_y(&t->mat, ind_vx(i, j), ind_vy(i, j), ind_sxy(i, j), ind_syy(i, j), &w);
		omega_y(&t->mat, ind_vx(i, j+1), ind_vy(i, j+1), ind_sxy(i, j+1), ind_syy(i, j+1), &w1);
		omega_y(&t->mat, ind_vx(i, j+2), ind_vy(i, j+2), ind_sxy(i, j+2), ind_syy(i, j+2), &w2);
		w_2 = w_1 = w;
	} else if (j == 1) {
		omega_y(&t->mat, ind_vx(i, j-1), ind_vy(i, j-1), ind_sxy(i, j-1), ind_syy(i, j-1), &w_1);
		omega_y(&t->mat, ind_vx(i, j), ind_vy(i, j), ind_sxy(i, j), ind_syy(i, j), &w);
		omega_y(&t->mat, ind_vx(i, j+1), ind_vy(i, j+1), ind_sxy(i, j+1), ind_syy(i, j+1), &w1);
		omega_y(&t->mat, ind_vx(i, j+2), ind_vy(i, j+2), ind_sxy(i, j+2), ind_syy(i, j+2), &w2);
		w_2 = w1;
	} else if (j == ny - 2) {
		omega_y(&t->mat, ind_vx(i, j-2), ind_vy(i, j-2), ind_sxy(i, j-2), ind_syy(i, j-2), &w_2);
		omega_y(&t->mat, ind_vx(i, j-1), ind_vy(i, j-1), ind_sxy(i, j-1), ind_syy(i, j-1), &w_1);
		omega_y(&t->mat, ind_vx(i, j), ind_vy(i, j), ind_sxy(i, j), ind_syy(i, j), &w);
		omega_y(&t->mat, ind_vx(i, j+1), ind_vy(i, j+1), ind_sxy(i, j+1), ind_syy(i, j+1), &w1);
		w2 = w1;
	} else if (j == ny - 1) {
		omega_y(&t->mat, ind_vx(i, j-2), ind_vy(i, j-2), ind_sxy(i, j-2), ind_syy(i, j-2), &w_2);
		omega_y(&t->mat, ind_vx(i, j-1), ind_vy(i, j-1), ind_sxy(i, j-1), ind_syy(i, j-1), &w_1);
		omega_y(&t->mat, ind_vx(i, j), ind_vy(i, j), ind_sxy(i, j), ind_syy(i, j), &w);
		w2 = w1 = w;
	}

	reconstruct(w_2, w_1, w, w1, w2, k1, k2, &d);
	inc_y(&t->mat, ind_vx(i,j), ind_vy(i,j), ind_sxx(i,j), ind_sxy(i,j), ind_syy(i,j), &vx, &vy, &sxx, &sxy, &syy, &d);
	grid = out_grid;
	write_new_values(&ind_vx(i,j), &ind_vy(i,j), &ind_sxx(i,j), &ind_sxy(i,j), &ind_syy(i,j), vx, vy, sxx, sxy, syy);
}

double le_step(le_task *task, int steps)
{
	// TODO place all this code in main

	int nx = task->n.x;
	int ny = task->n.y;

	// set sizes of blocks on gpu
	int threads_width = 32;
	dim3 threadsPerBlock(threads_width, threads_width);
    dim3 blocksPerGrid((nx + threads_width - 1) / threads_width, (ny  + threads_width - 1) / threads_width);

    int grid_size = sizeof(real) * NODE_SIZE * nx * ny;
	int task_size = sizeof(le_task);
	le_task* d_task;
	real* d_grid1;
	real* d_grid2;

	double t;
    
    // allocate memory on gpu
	gpuErrchk(cudaMalloc(&d_task, task_size));
	gpuErrchk(cudaMalloc(&d_grid1, grid_size));
	gpuErrchk(cudaMalloc(&d_grid2, grid_size));
	gpuErrchk(cudaMemcpy(d_task, task, task_size, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_grid1, task->grid, grid_size, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_grid2, task->grid, grid_size, cudaMemcpyHostToDevice));

	cudaDeviceSynchronize();
	t = timer();
	// run kernel
	for (int i = 0; i < steps; i++) {
		le_step_x<<<blocksPerGrid, threadsPerBlock>>>(d_task, d_grid1, d_grid2);
		gpuErrchk( cudaPeekAtLastError() );
		le_step_y<<<blocksPerGrid, threadsPerBlock>>>(d_task, d_grid2, d_grid1);
		gpuErrchk( cudaPeekAtLastError() );
	}
	gpuErrchk(cudaDeviceSynchronize());
	t = timer() - t;

	// drop data to host
	gpuErrchk(cudaMemcpy(task->grid, d_grid1, grid_size, cudaMemcpyDeviceToHost));
	cudaFree(d_grid1);
	cudaFree(d_grid2);
	cudaFree(d_task);

	return t;
}
