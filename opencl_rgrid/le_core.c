#include <malloc.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#define USE_OPENCL

#include "CL/cl.h"

#include "rgrid/darraycontainer.h"
#include "rgrid/clwrapper.h"
#include "rgrid/clutils.h"

#include "le_types.h"
#include "le_core.h"

#define MAX_ENTRIES 10
#define BUF_SIZE 160000 // for log buffer
#define SOURCE_BUF_SIZE 16000

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

#define vnorm(v) (sqrt(v.x * v.x + v.y * v.y))

inline real le_min(real a, real b) { return a > b ? b : a; }
inline real le_max(real a, real b) { return a > b ? a : b; }
inline real le_max3(real a, real b, real c) { return le_max(a, le_max(b, c)); }

#ifdef USE_DOUBLE

#define TVD2_EPS 1e-6

#define limiter_minmod(r) (le_max(0.0, le_min(1.0, (r))))
#define limiter_cir(r) (0.0)
#define limiter_superbee(r) (le_max3(0.0, le_min(1.0, 2.0 * r), le_min(2.0, r)))

const char* options = "-DUSE_DOUBLE";

#else

#define TVD2_EPS 1e-6f

#define limiter_minmod(r) (le_max(0.0f, le_min(1.0f, (r))))
#define limiter_cir(r) (0.0f)
#define limiter_superbee(r) (le_max3(0.0f, le_min(1.0f, 2.0f * r), le_min(2.0f, r)))

const char* options = "";

#endif

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

inline void inc_x(const le_material *m, const real vx, const real vy, const real sxx, const real sxy, const real syy,
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

inline void inc_y(const le_material *m, const real vx, const real vy, const real sxx, const real sxy, const real syy,
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

inline void reconstruct(const le_w ppu, const le_w pu, const le_w u, const le_w nu, const le_w nnu, const real k1, const real k2, le_w *d)
{
	d->w1 = tvd2(k1, ppu.w1, pu.w1, u.w1, nu.w1); // c1
	d->w2 = tvd2(k1, nnu.w2, nu.w2, u.w2, pu.w2); // -c1
	d->w3 = tvd2(k2, ppu.w3, pu.w3, u.w3, nu.w3); // c2
	d->w4 = tvd2(k2, nnu.w4, nu.w4, u.w4, pu.w4); // -c2
}

inline real g_ind(const real* grid, int i, int j, const int nx, const int ny, const int node)
{
	// TODO it works only with SOA

	if (i < 0) i = 0;
	if (i >= nx) i = nx - 1;
	if (j < 0) j = 0;
	if (j >= ny) j = ny - 1;
	return (*(grid + (i) + (j) * nx + node * nx * ny));
}

void loadFromFile(const char* filename, char* source, int max_size)
{
	FILE* input = fopen(filename, "r");
	if (input == NULL) {
		printf("Can't open %s\n", filename);
		exit(-1);
	}
	int pos = 0;
	while (!feof(input)) {
		source[pos++] = fgetc(input);
		if (pos == max_size) { 
			printf("Can't load source from %s, file so big\n", filename);
			exit(-1);
		}
	}
	source[pos-1] = '\0';
	fclose(input);
}

double le_step(le_task *task, int steps)
{
	const int nx = task->n.x;
	const int ny = task->n.y;
	
	const int grid_size = sizeof(real) * NODE_SIZE * nx * ny;
	
	rgrid::DArray<real, int> da;
	da.resize(nx, ny, 1, 
	         nx, ny, 1,
	         0, 0, 0,
	         2, 2, 0,
	         NODE_SIZE);
	
	for (int cn = 0; cn != NODE_SIZE; ++cn)
	for (int j = 0; j != ny; ++j)
	for (int i = 0; i != nx; ++i) {
		//da[(i + 2) + (j + 2) * (nx + 4) + cn * (nx + 4) * (ny + 4)] = task->grid[i + j * nx + cn * nx * ny];
		da(i, j, cn) = task->grid[i + j * nx + cn * nx * ny];
	}
	
	const clwrapper::CLWrapper& clw = clwrapper::CLWrapper::instance();
	
	assert(clw.getPlatformsNum() == 1);
	assert(clw.getDevicesNum() != 0);
	
	rgrid::DArrayContainer<real, int> dac(da, 1, clw.getDevicesNum(), 1);
	
	const real k1x = task->dt * task->mat.c1 / task->h.x;
	const real k2x = task->dt * task->mat.c2 / task->h.x;
	const real k1y = task->dt * task->mat.c1 / task->h.y;
	const real k2y = task->dt * task->mat.c2 / task->h.y;

	// set sizes of blocks on gpu
	size_t block_dim_stepx_x = 256;
	size_t block_dim_stepx_y = 1;
	size_t local_work_size_x[2] = {
		block_dim_stepx_x, 
		block_dim_stepx_y};
	size_t global_work_size_x[2] = {
		(nx + block_dim_stepx_x - 1) / block_dim_stepx_x * block_dim_stepx_x,
		(ny  + block_dim_stepx_y - 1) / block_dim_stepx_y * block_dim_stepx_y};
	size_t sharedMemSizeX = (block_dim_stepx_x + 4) * block_dim_stepx_y * NODE_SIZE * sizeof(real);
	
	size_t block_dim_stepy_x = 16;
	size_t block_dim_stepy_y = 16;
	size_t local_work_size_y[2] = {block_dim_stepy_x, block_dim_stepy_y};
	size_t global_work_size_y[2] = {(nx + block_dim_stepy_x - 1) / block_dim_stepy_x * block_dim_stepy_x, 
		(ny  + block_dim_stepy_y - 1) / block_dim_stepy_y * block_dim_stepy_y};
	size_t sharedMemSizeY = block_dim_stepy_x * (block_dim_stepy_y + 4) * NODE_SIZE * sizeof(real);
	
	
	cl_int err;
	// load program source from file
	char* source = (char*) malloc(sizeof(char) * SOURCE_BUF_SIZE);
	assert(source);
	loadFromFile("kernel.cl", source, SOURCE_BUF_SIZE);
	
	// create program
	cl_program program = clCreateProgramWithSource(clw.getContext(), 1, (const char **) &source, NULL, &err);
	CHECK_CL_ERROR(err);
	
	// build program
	err = clBuildProgram(program, 0, NULL, options, NULL, NULL);
	
	// print build info
	char buf[BUF_SIZE];
	cl_build_status build_status;
	CHECK_CL_ERROR(clGetProgramBuildInfo(program, clw.getDevice(0), CL_PROGRAM_BUILD_LOG, BUF_SIZE * sizeof(char), buf, NULL));
	CHECK_CL_ERROR(clGetProgramBuildInfo(program, clw.getDevice(0), CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), &build_status, NULL));
	printf("build log:\n%s\n", buf);
	if (CL_BUILD_SUCCESS == build_status) printf("Build program: OK\n");
	else {
		printf("Build program: FAIL\n");
		exit(-1);
	}

	// check clBuildProgram
	CHECK_CL_ERROR(err);
	
	cl_kernel kernel_step_x[MAX_ENTRIES], kernel_step_y[MAX_ENTRIES];
	
	for (int dev = 0; dev != clw.getDevicesNum(); ++dev) {
		rgrid::DArray<real, int>& dap = dac.getDArrayPart(dev);
		dap.setCLContext(clw.getContext());
		dap.setCLCQ(clw.getCommandQueue(dev));
		dap.clHtoD();
		
		// create kernel
		kernel_step_x[dev] = clCreateKernel(program, "le_step_x", &err);
		CHECK_CL_ERROR(err);
		kernel_step_y[dev] = clCreateKernel(program, "le_step_y", &err);
		CHECK_CL_ERROR(err);
		
		int lnx = dap.localSize(rgrid::X);
		int lny = dap.localSize(rgrid::Y);
		
		// set kernel args
		CHECK_CL_ERROR(clSetKernelArg(kernel_step_x[dev], 0, sizeof(lnx), &lnx));
		CHECK_CL_ERROR(clSetKernelArg(kernel_step_x[dev], 1, sizeof(lny), &lny));
		CHECK_CL_ERROR(clSetKernelArg(kernel_step_x[dev], 2, sizeof(k1x), &k1x));
		CHECK_CL_ERROR(clSetKernelArg(kernel_step_x[dev], 3, sizeof(k2x), &k2x));
		CHECK_CL_ERROR(clSetKernelArg(kernel_step_x[dev], 4, sizeof(task->mat), &(task->mat)));
		CHECK_CL_ERROR(clSetKernelArg(kernel_step_x[dev], 7, sharedMemSizeX, NULL));
		
		CHECK_CL_ERROR(clSetKernelArg(kernel_step_y[dev], 0, sizeof(lnx), &lnx));
		CHECK_CL_ERROR(clSetKernelArg(kernel_step_y[dev], 1, sizeof(lny), &lny));
		CHECK_CL_ERROR(clSetKernelArg(kernel_step_y[dev], 2, sizeof(k1y), &k1y));
		CHECK_CL_ERROR(clSetKernelArg(kernel_step_y[dev], 3, sizeof(k2y), &k2y));
		CHECK_CL_ERROR(clSetKernelArg(kernel_step_y[dev], 4, sizeof(task->mat), &(task->mat)));
		CHECK_CL_ERROR(clSetKernelArg(kernel_step_y[dev], 7, sharedMemSizeY, NULL));
	}
	
	double t = timer();
	// add kernel execution in queue
	for (int i = 0; i < steps; i++)
	{
		dac.syncCL();
		dac.fillGhostCL();
		for (int dev = 0; dev != clw.getDevicesNum(); ++dev) {
			rgrid::DArray<real, int>& dap = dac.getDArrayPart(dev);
			cl_mem d_grid1 = dap.getCLBuffer();
			cl_mem d_grid2 = dap.getCLBuffer2();
			CHECK_CL_ERROR(clSetKernelArg(kernel_step_x[dev], 5, sizeof(d_grid1), &d_grid1));
			CHECK_CL_ERROR(clSetKernelArg(kernel_step_x[dev], 6, sizeof(d_grid2), &d_grid2));
			CHECK_CL_ERROR(clEnqueueNDRangeKernel(clw.getCommandQueue(dev), kernel_step_x[dev], 2, NULL, global_work_size_x, local_work_size_x, 0, NULL, NULL));
			dap.swapCLBuffers();
		}
		dac.syncCL();
		dac.fillGhostCL();
		for (int dev = 0; dev != clw.getDevicesNum(); ++dev) {
			rgrid::DArray<real, int>& dap = dac.getDArrayPart(dev);
			cl_mem d_grid1 = dap.getCLBuffer();
			cl_mem d_grid2 = dap.getCLBuffer2();
			CHECK_CL_ERROR(clSetKernelArg(kernel_step_y[dev], 5, sizeof(d_grid1), &d_grid1));
			CHECK_CL_ERROR(clSetKernelArg(kernel_step_y[dev], 6, sizeof(d_grid2), &d_grid2));
			CHECK_CL_ERROR(clEnqueueNDRangeKernel(clw.getCommandQueue(dev), kernel_step_y[dev], 2, NULL, global_work_size_y, local_work_size_y, 0, NULL, NULL));
			dap.swapCLBuffers();
		}
	}
	for (int dev = 0; dev != clw.getDevicesNum(); ++dev) {
		CHECK_CL_ERROR(clFinish(clw.getCommandQueue(dev)));
	}
	t = timer() - t;

	for (int dev = 0; dev != clw.getDevicesNum(); ++dev) {
		rgrid::DArray<real, int>& dap = dac.getDArrayPart(dev);
		// return data to host
		dap.clDtoH();
		// release kernel
		CHECK_CL_ERROR(clReleaseKernel(kernel_step_x[dev]));
		CHECK_CL_ERROR(clReleaseKernel(kernel_step_y[dev]));
	}
	
	// release program
	CHECK_CL_ERROR(clReleaseProgram(program));
	
	free(source);
	
	dac.getDArray(da);
	for (int cn = 0; cn != NODE_SIZE; ++cn)
	for (int j = 0; j != ny; ++j)
	for (int i = 0; i != nx; ++i) {
		//task->grid[i + j * nx + cn * nx * ny] = da[(i + 2) + (j + 2) * (nx + 4) + cn * (nx + 4) * (ny + 4)];
		task->grid[i + j * nx + cn * nx * ny] = da(i, j, cn);
	}
	
	return t;
}


