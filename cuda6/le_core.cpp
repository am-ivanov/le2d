#include <malloc.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include <cuda.h>
#include <cuda_runtime.h>

#include "le_core.h"

#define s_indx_vx(i, j) (*(shared_grid + ((i)+2) + (j) * (blockDim.x + 4) + 0 * (blockDim.x + 4) * blockDim.y))
#define s_indx_vy(i, j) (*(shared_grid + ((i)+2) + (j) * (blockDim.x + 4) + 1 * (blockDim.x + 4) * blockDim.y))
#define s_indx_sxx(i, j) (*(shared_grid + ((i)+2) + (j) * (blockDim.x + 4) + 2 * (blockDim.x + 4) * blockDim.y))
#define s_indx_sxy(i, j) (*(shared_grid + ((i)+2) + (j) * (blockDim.x + 4) + 3 * (blockDim.x + 4) * blockDim.y))
#define s_indx_syy(i, j) (*(shared_grid + ((i)+2) + (j) * (blockDim.x + 4) + 4 * (blockDim.x + 4) * blockDim.y))

#define s_indy_vx(i, j) (*(shared_grid + (i) + ((j)+2) * blockDim.x + 0 * blockDim.x * (blockDim.y + 4)))
#define s_indy_vy(i, j) (*(shared_grid + (i) + ((j)+2) * blockDim.x + 1 * blockDim.x * (blockDim.y + 4)))
#define s_indy_sxx(i, j) (*(shared_grid + (i) + ((j)+2) * blockDim.x + 2 * blockDim.x * (blockDim.y + 4)))
#define s_indy_sxy(i, j) (*(shared_grid + (i) + ((j)+2) * blockDim.x + 3 * blockDim.x * (blockDim.y + 4)))
#define s_indy_syy(i, j) (*(shared_grid + (i) + ((j)+2) * blockDim.x + 4 * blockDim.x * (blockDim.y + 4)))

#define MAX_DEVICE_NUM 10

#ifdef USE_DOUBLE

#define TVD2_EPS 1e-6
#define CONST_ZERO 0.0
#define CONST_HALF 0.5
#define CONST_ONE 1.0
#define CONST_TWO 2.0

#else // FLOAT

#define TVD2_EPS 1e-6f
#define CONST_ZERO 0.0f
#define CONST_HALF 0.5f
#define CONST_ONE 1.0f
#define CONST_TWO 2.0f

#endif

#define vnorm(v) (sqrt(v.x * v.x + v.y * v.y))

__device__ __inline__ real le_min(real a, real b) { return a > b ? b : a; }
__device__ __inline__ real le_max(real a, real b) { return a > b ? a : b; }
__device__ __inline__ real le_max3(real a, real b, real c) { return le_max(a, le_max(b, c)); }

#define limiter_minmod(r) (le_max(CONST_ZERO, le_min(CONST_ONE, (r))))
#define limiter_cir(r) (CONST_ZERO)
#define limiter_superbee(r) (le_max3(CONST_ZERO, le_min(CONST_ONE, CONST_TWO * r), le_min(CONST_TWO, r)))

#define limiter limiter_superbee

__device__ __inline__ real tvd2(const real c, const real u_2, const real u_1, const real u, const real u1)
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

	const real k = CONST_HALF * (CONST_ONE - c);

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

__device__ __inline__ void omega_x(const le_material *m, const real vx, const real vy, const real sxx, const real sxy, le_w *w)
{
	const real N00T = sxx * m->irhoc1;
	const real N01T = sxy * m->irhoc2;

	w->w1 = vx - N00T;
	w->w2 = vx + N00T;
	w->w3 = vy - N01T;
	w->w4 = vy + N01T;
}

__device__ __inline__ void omega_y(const le_material *m, const real vx, const real vy, const real sxy, const real syy, le_w *w)
{
	const real N00T = syy * m->irhoc1;
	const real N01T = sxy * m->irhoc2;

	w->w1 = vy - N00T;
	w->w2 = vy + N00T;
	w->w3 = vx - N01T;
	w->w4 = vx + N01T;
}

__device__ __inline__ void inc_x(const le_material *m, const real vx, const real vy, const real sxx, const real sxy, const real syy,
		real* new_vx, real* new_vy, real* new_sxx, real* new_sxy, real* new_syy, const le_w *d)
{
	const real d1 = CONST_HALF * d->w1;
	const real d2 = CONST_HALF * d->w2;
	const real d3 = CONST_HALF * d->w3;
	const real d4 = CONST_HALF * d->w4;

	*new_vx = vx + d1 + d2;
	*new_vy = vy + d3 + d4;

	*new_syy = syy + (d2 - d1) * m->rhoc3;
	*new_sxx = sxx + (d2 - d1) * m->rhoc1;
	*new_sxy = sxy + (d4 - d3) * m->rhoc2;
}

__device__ __inline__ void inc_y(const le_material *m, const real vx, const real vy, const real sxx, const real sxy, const real syy,
		real* new_vx, real* new_vy, real* new_sxx, real* new_sxy, real* new_syy, const le_w *d)
{
	const real d1 = CONST_HALF * d->w1;
	const real d2 = CONST_HALF * d->w2;
	const real d3 = CONST_HALF * d->w3;
	const real d4 = CONST_HALF * d->w4;

	*new_vy = vy + d1 + d2;
	*new_vx = vx + d3 + d4;

	*new_syy = syy + (d2 - d1) * m->rhoc1;
	*new_sxx = sxx + (d2 - d1) * m->rhoc3;
	*new_sxy = sxy + (d4 - d3) * m->rhoc2;
}

__device__ __inline__ void reconstruct(const le_w ppu, const le_w pu, const le_w u, const le_w nu, const le_w nnu, const real k1, const real k2, le_w *d)
{
	d->w1 = tvd2(k1, ppu.w1, pu.w1, u.w1, nu.w1); // c1
	d->w2 = tvd2(k1, nnu.w2, nu.w2, u.w2, pu.w2); // -c1
	d->w3 = tvd2(k2, ppu.w3, pu.w3, u.w3, nu.w3); // c2
	d->w4 = tvd2(k2, nnu.w4, nu.w4, u.w4, pu.w4); // -c2
}

__device__ __inline__ real 
g_ind_x(const real* grid, int i, int j, const int nx, const int ny, const int node)
{
	if (i < 0) i = 0;
	if (i >= nx) i = nx - 1;
	if (j < 0) j = 0;
	if (j >= ny) j = ny - 1;
	return (*(grid + (i) + (j) * nx + node * nx * ny));
}

__device__ __inline__ real 
g_ind_y(const real* grid, int i, int j, const int nx, const int ny, const int node, const real* grid_top, const real* grid_bot)
{
	if (i < 0) i = 0;
	if (i >= nx) i = nx - 1;
	if (j < -2) j = 0;
	else if (j < 0) {
		// j == -2 or -1
		return (*(grid_bot + (i) + (j + 2) * nx + node * nx * 2));
	} 
	if (j >= ny + 2) j = ny - 1;
	else if (j >= ny) {
		// j == (ny) or (ny + 1)
		return (*(grid_top + (i) + (j - ny) * nx + node * nx * 2));
	}
	return (*(grid + (i) + (j) * nx + node * nx * ny));
}

extern __shared__ real shared_grid[];

__global__ void le_step_x(const int nx, const int ny, const real k1, const real k2, const le_material mat,
		const real* __restrict__ in_grid, real* __restrict__ grid)
{
	const int i = threadIdx.x + blockDim.x * blockIdx.x;
	const int j = threadIdx.y + blockDim.y * blockIdx.y;
	const int li = threadIdx.x;
	const int lj = threadIdx.y;

	real vx, vy, sxx, sxy, syy;
	le_w w_2, w_1, w, w1, w2, d;

	s_indx_vx(li, lj) = g_ind_x(in_grid, i, j, nx, ny, 0);
	s_indx_vy(li, lj) = g_ind_x(in_grid, i, j, nx, ny, 1);
	s_indx_sxx(li, lj) = g_ind_x(in_grid, i, j, nx, ny, 2);
	s_indx_sxy(li, lj) = g_ind_x(in_grid, i, j, nx, ny, 3);
	s_indx_syy(li, lj) = g_ind_x(in_grid, i, j, nx, ny, 4);
	if (li < 2) {
		s_indx_vx(li - 2, lj) = g_ind_x(in_grid, i - 2, j, nx, ny, 0);
		s_indx_vy(li - 2, lj) = g_ind_x(in_grid, i - 2, j, nx, ny, 1);
		s_indx_sxx(li - 2, lj) = g_ind_x(in_grid, i - 2, j, nx, ny, 2);
		s_indx_sxy(li - 2, lj) = g_ind_x(in_grid, i - 2, j, nx, ny, 3);
		s_indx_syy(li - 2, lj) = g_ind_x(in_grid, i - 2, j, nx, ny, 4);
	} else if (li >= blockDim.x - 2) {
		s_indx_vx(li + 2, lj) = g_ind_x(in_grid, i + 2, j, nx, ny, 0);
		s_indx_vy(li + 2, lj) = g_ind_x(in_grid, i + 2, j, nx, ny, 1);
		s_indx_sxx(li + 2, lj) = g_ind_x(in_grid, i + 2, j, nx, ny, 2);
		s_indx_sxy(li + 2, lj) = g_ind_x(in_grid, i + 2, j, nx, ny, 3);
		s_indx_syy(li + 2, lj) = g_ind_x(in_grid, i + 2, j, nx, ny, 4);
	}

	__syncthreads();

	if (i >= nx || j >= ny) return;

	omega_x(&mat, s_indx_vx(li-2, lj), s_indx_vy(li-2, lj), s_indx_sxx(li-2, lj), s_indx_sxy(li-2, lj), &w_2);
	omega_x(&mat, s_indx_vx(li-1, lj), s_indx_vy(li-1, lj), s_indx_sxx(li-1, lj), s_indx_sxy(li-1, lj), &w_1);
	omega_x(&mat, s_indx_vx(li, lj), s_indx_vy(li, lj), s_indx_sxx(li, lj), s_indx_sxy(li, lj), &w);
	omega_x(&mat, s_indx_vx(li+1, lj), s_indx_vy(li+1, lj), s_indx_sxx(li+1, lj), s_indx_sxy(li+1, lj), &w1);
	omega_x(&mat, s_indx_vx(li+2, lj), s_indx_vy(li+2, lj), s_indx_sxx(li+2, lj), s_indx_sxy(li+2, lj), &w2);

	reconstruct(w_2, w_1, w, w1, w2, k1, k2, &d);
	inc_x(&mat, s_indx_vx(li,lj), s_indx_vy(li,lj), s_indx_sxx(li,lj), s_indx_sxy(li,lj), s_indx_syy(li,lj), &vx, &vy, &sxx, &sxy, &syy, &d);

	ind_vx(i,j) = vx;
	ind_vy(i,j) = vy;
	ind_sxx(i,j) = sxx;
	ind_sxy(i,j) = sxy;
	ind_syy(i,j) = syy;
}

__global__ void le_step_y(const int nx, const int ny, const real k1, const real k2, const le_material mat,
		const real* __restrict__ in_grid, real* __restrict__ grid, const real* __restrict__ grid_top, const real* __restrict__ grid_bot)
{
	const int i = threadIdx.x + blockDim.x * blockIdx.x;
	const int j = threadIdx.y + blockDim.y * blockIdx.y;
	const int li = threadIdx.x;
	const int lj = threadIdx.y;

	real vx, vy, sxx, sxy, syy;
	le_w w_2, w_1, w, w1, w2, d;

	s_indy_vx(li, lj) = g_ind_y(in_grid, i, j, nx, ny, 0, grid_top, grid_bot);
	s_indy_vy(li, lj) = g_ind_y(in_grid, i, j, nx, ny, 1, grid_top, grid_bot);
	s_indy_sxx(li, lj) = g_ind_y(in_grid, i, j, nx, ny, 2, grid_top, grid_bot);
	s_indy_sxy(li, lj) = g_ind_y(in_grid, i, j, nx, ny, 3, grid_top, grid_bot);
	s_indy_syy(li, lj) = g_ind_y(in_grid, i, j, nx, ny, 4, grid_top, grid_bot);
	if (lj < 2) {
		s_indy_vx(li, lj - 2) = g_ind_y(in_grid, i, j - 2, nx, ny, 0, grid_top, grid_bot);
		s_indy_vy(li, lj - 2) = g_ind_y(in_grid, i, j - 2, nx, ny, 1, grid_top, grid_bot);
		s_indy_sxx(li, lj - 2) = g_ind_y(in_grid, i, j - 2, nx, ny, 2, grid_top, grid_bot);
		s_indy_sxy(li, lj - 2) = g_ind_y(in_grid, i, j - 2, nx, ny, 3, grid_top, grid_bot);
		s_indy_syy(li, lj - 2) = g_ind_y(in_grid, i, j - 2, nx, ny, 4, grid_top, grid_bot);
	} else if (lj >= blockDim.y - 2) {
		s_indy_vx(li, lj + 2) = g_ind_y(in_grid, i, j + 2, nx, ny, 0, grid_top, grid_bot);
		s_indy_vy(li, lj + 2) = g_ind_y(in_grid, i, j + 2, nx, ny, 1, grid_top, grid_bot);
		s_indy_sxx(li, lj + 2) = g_ind_y(in_grid, i, j + 2, nx, ny, 2, grid_top, grid_bot);
		s_indy_sxy(li, lj + 2) = g_ind_y(in_grid, i, j + 2, nx, ny, 3, grid_top, grid_bot);
		s_indy_syy(li, lj + 2) = g_ind_y(in_grid, i, j + 2, nx, ny, 4, grid_top, grid_bot);
	}

	__syncthreads();

	if (i >= nx || j >= ny) return;

	omega_y(&mat, s_indy_vx(li, lj-2), s_indy_vy(li, lj-2), s_indy_sxy(li, lj-2), s_indy_syy(li, lj-2), &w_2);
	omega_y(&mat, s_indy_vx(li, lj-1), s_indy_vy(li, lj-1), s_indy_sxy(li, lj-1), s_indy_syy(li, lj-1), &w_1);
	omega_y(&mat, s_indy_vx(li, lj), s_indy_vy(li, lj), s_indy_sxy(li, lj), s_indy_syy(li, lj), &w);
	omega_y(&mat, s_indy_vx(li, lj+1), s_indy_vy(li, lj+1), s_indy_sxy(li, lj+1), s_indy_syy(li, lj+1), &w1);
	omega_y(&mat, s_indy_vx(li, lj+2), s_indy_vy(li, lj+2), s_indy_sxy(li, lj+2), s_indy_syy(li, lj+2), &w2);

	reconstruct(w_2, w_1, w, w1, w2, k1, k2, &d);
	inc_y(&mat, s_indy_vx(li,lj), s_indy_vy(li,lj), s_indy_sxx(li,lj), s_indy_sxy(li,lj), s_indy_syy(li,lj), &vx, &vy, &sxx, &sxy, &syy, &d);

	ind_vx(i,j) = vx;
	ind_vy(i,j) = vy;
	ind_sxx(i,j) = sxx;
	ind_sxy(i,j) = sxy;
	ind_syy(i,j) = syy;
}

double le_step(le_task *task, int steps)
{
	const int nx = task->n.x;
	const int ny = task->n.y;

	const real k1x = task->dt * task->mat.c1 / task->h.x;
	const real k2x = task->dt * task->mat.c2 / task->h.x;
	const real k1y = task->dt * task->mat.c1 / task->h.y;
	const real k2y = task->dt * task->mat.c2 / task->h.y;

	// get num of devices
	int deviceCount;
	checkCudaErrors(cudaGetDeviceCount(&deviceCount));
	
	// calcualte grid size	
	int dev_ny = ny / deviceCount;
	if (ny % deviceCount != 0) {
		printf("wrong ny size:\n");
		printf("ny %% deviceCount == %d\n", ny % deviceCount);
		exit(-1);
	}
	real* d_grid1[MAX_DEVICE_NUM];
	real* d_grid2[MAX_DEVICE_NUM];
	real* d_top[MAX_DEVICE_NUM];
	real* d_bot[MAX_DEVICE_NUM];
	
	// set sizes of blocks on gpu
	int block_dim_stepx_x = 256;
	int block_dim_stepx_y = 1;
	dim3 threadsPerBlockX(block_dim_stepx_x, block_dim_stepx_y);
	dim3 blocksPerGridX((nx + block_dim_stepx_x - 1) / block_dim_stepx_x, (dev_ny  + block_dim_stepx_y - 1) / block_dim_stepx_y);
	int sharedMemSizeX = (block_dim_stepx_x + 4) * block_dim_stepx_y * NODE_SIZE * sizeof(real);
	
	int block_dim_stepy_x = 16;
	int block_dim_stepy_y = 16;
	dim3 threadsPerBlockY(block_dim_stepy_x, block_dim_stepy_y);
	dim3 blocksPerGridY((nx + block_dim_stepy_x - 1) / block_dim_stepy_x, (dev_ny  + block_dim_stepy_y - 1) / block_dim_stepy_y);
	int sharedMemSizeY = block_dim_stepy_x * (block_dim_stepy_y + 4) * NODE_SIZE * sizeof(real);

	// enable peer access from all to all devices
	for (int i = 0; i < deviceCount; i++) {
		checkCudaErrors(cudaSetDevice(i));
		for (int j = 0; j < deviceCount; j++) {
			if (j == i) continue;
			int canAccessPeer;
			checkCudaErrors(cudaDeviceCanAccessPeer(&canAccessPeer, i, j));
			if (canAccessPeer == 0) {
				printf("Peer acces from %d to %d is not allowed\n", i, j);
			}
			else 
				checkCudaErrors(cudaDeviceEnablePeerAccess(j, 0));
		}
	}
	
	// synchronize all devices
	for (int i = 0; i < deviceCount; i++) {
		checkCudaErrors(cudaSetDevice(i));
		checkCudaErrors(cudaDeviceSynchronize());
	}
	
	double t = timer();
	
	// allocate memory on gpu
	for (int i = 0; i < deviceCount; i++)
	{
		checkCudaErrors(cudaSetDevice(i));
		checkCudaErrors(cudaMalloc(&(d_grid1[i]), sizeof(real) * NODE_SIZE * nx * dev_ny));
		checkCudaErrors(cudaMalloc(&(d_grid2[i]), sizeof(real) * NODE_SIZE * nx * dev_ny));
		checkCudaErrors(cudaMalloc(&(d_top[i]), sizeof(real) * NODE_SIZE * nx * 2));
		checkCudaErrors(cudaMalloc(&(d_bot[i]), sizeof(real) * NODE_SIZE * nx * 2));
		for (int j = 0; j < NODE_SIZE; j++) {
			checkCudaErrors(cudaMemcpy(d_grid1[i] + j * nx * dev_ny, task->grid + j * nx * ny + nx * dev_ny * i, sizeof(real) * nx * dev_ny, cudaMemcpyHostToDevice));
		}
	}
	
	for (int i = 0; i < steps; i++) {
		for (int j = 0; j < deviceCount; j++) {
			checkCudaErrors(cudaSetDevice(j));
			le_step_x<<<blocksPerGridX, threadsPerBlockX, sharedMemSizeX>>>(nx, dev_ny, k1x, k2x, task->mat, d_grid1[j], d_grid2[j]);
			checkCudaErrors(cudaPeekAtLastError());
		}
		// copy memory between devices before step y
		for (int j = 0; j < deviceCount; j++)
		{
			checkCudaErrors(cudaSetDevice(j));
			if (j == 0) {
				for (int k = 0; k < NODE_SIZE; k++) {
					checkCudaErrors(cudaMemcpy(d_bot[j] + k * nx * 2, d_grid2[j] + k * nx * dev_ny, sizeof(real) * nx, cudaMemcpyDeviceToDevice));
					checkCudaErrors(cudaMemcpy(d_bot[j] + k * nx * 2 + nx, d_grid2[j] + k * nx * dev_ny, sizeof(real) * nx, cudaMemcpyDeviceToDevice));
				}
			} else {
				for (int k = 0; k < NODE_SIZE; k++) {
					checkCudaErrors(cudaMemcpyPeer(d_bot[j] + k * nx * 2, j, d_grid2[j - 1] + k * nx * dev_ny + nx * (dev_ny - 2), j - 1, sizeof(real) * nx * 2));
				}
			}
			if (j == deviceCount - 1) {
				for (int k = 0; k < NODE_SIZE; k++) {
					checkCudaErrors(cudaMemcpy(d_top[j] + k * nx * 2, d_grid2[j] + k * nx * dev_ny + nx * (dev_ny - 1), sizeof(real) * nx, cudaMemcpyDeviceToDevice));
					checkCudaErrors(cudaMemcpy(d_top[j] + k * nx * 2 + nx, d_grid2[j] + k * nx * dev_ny + nx * (dev_ny - 1), sizeof(real) * nx, cudaMemcpyDeviceToDevice));
				}
			} else {
				for (int k = 0; k < NODE_SIZE; k++) {
					checkCudaErrors(cudaMemcpyPeer(d_top[j] + k * nx * 2, j, d_grid2[j + 1] + k * nx * dev_ny, j + 1, sizeof(real) * nx * 2));
				}
			}
		}
		for (int j = 0; j < deviceCount; j++) {
			checkCudaErrors(cudaSetDevice(j));
			le_step_y<<<blocksPerGridY, threadsPerBlockY, sharedMemSizeY>>>(nx, dev_ny, k1y, k2y, task->mat, d_grid2[j], d_grid1[j], d_top[j], d_bot[j]);
			checkCudaErrors(cudaPeekAtLastError());
		}
	}

	// return data to host
	for (int i = 0; i < deviceCount; i++)
	{
		checkCudaErrors(cudaSetDevice(i));
		for (int j = 0; j < NODE_SIZE; j++) {
			checkCudaErrors(cudaMemcpy(task->grid + j * nx * ny + nx * dev_ny * i, d_grid1[i] + j * nx * dev_ny, sizeof(real) * nx * dev_ny, cudaMemcpyDeviceToHost));
		}
		checkCudaErrors(cudaFree(d_grid1[i]));
		checkCudaErrors(cudaFree(d_grid2[i]));
		checkCudaErrors(cudaFree(d_top[i]));
		checkCudaErrors(cudaFree(d_bot[i]));
	}
	
	// synchronize all devices
	for (int i = 0; i < deviceCount; i++) {
		checkCudaErrors(cudaSetDevice(i));
		checkCudaErrors(cudaDeviceSynchronize());
	}
	
	t = timer() - t;

	return t;
}



