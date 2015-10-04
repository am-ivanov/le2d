#include <malloc.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "le_types.h"
#include "le_core.h"

#include "CL/opencl.h"

const char *cl_error_to_str(cl_int e)
{
  switch (e)
  {
    case CL_SUCCESS: return "success";
    case CL_DEVICE_NOT_FOUND: return "device not found";
    case CL_DEVICE_NOT_AVAILABLE: return "device not available";
#if !(defined(CL_PLATFORM_NVIDIA) && CL_PLATFORM_NVIDIA == 0x3001)
    case CL_COMPILER_NOT_AVAILABLE: return "device compiler not available";
#endif
    case CL_MEM_OBJECT_ALLOCATION_FAILURE: return "mem object allocation failure";
    case CL_OUT_OF_RESOURCES: return "out of resources";
    case CL_OUT_OF_HOST_MEMORY: return "out of host memory";
    case CL_PROFILING_INFO_NOT_AVAILABLE: return "profiling info not available";
    case CL_MEM_COPY_OVERLAP: return "mem copy overlap";
    case CL_IMAGE_FORMAT_MISMATCH: return "image format mismatch";
    case CL_IMAGE_FORMAT_NOT_SUPPORTED: return "image format not supported";
    case CL_BUILD_PROGRAM_FAILURE: return "build program failure";
    case CL_MAP_FAILURE: return "map failure";

    case CL_INVALID_VALUE: return "invalid value";
    case CL_INVALID_DEVICE_TYPE: return "invalid device type";
    case CL_INVALID_PLATFORM: return "invalid platform";
    case CL_INVALID_DEVICE: return "invalid device";
    case CL_INVALID_CONTEXT: return "invalid context";
    case CL_INVALID_QUEUE_PROPERTIES: return "invalid queue properties";
    case CL_INVALID_COMMAND_QUEUE: return "invalid command queue";
    case CL_INVALID_HOST_PTR: return "invalid host ptr";
    case CL_INVALID_MEM_OBJECT: return "invalid mem object";
    case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR: return "invalid image format descriptor";
    case CL_INVALID_IMAGE_SIZE: return "invalid image size";
    case CL_INVALID_SAMPLER: return "invalid sampler";
    case CL_INVALID_BINARY: return "invalid binary";
    case CL_INVALID_BUILD_OPTIONS: return "invalid build options";
    case CL_INVALID_PROGRAM: return "invalid program";
    case CL_INVALID_PROGRAM_EXECUTABLE: return "invalid program executable";
    case CL_INVALID_KERNEL_NAME: return "invalid kernel name";
    case CL_INVALID_KERNEL_DEFINITION: return "invalid kernel definition";
    case CL_INVALID_KERNEL: return "invalid kernel";
    case CL_INVALID_ARG_INDEX: return "invalid arg index";
    case CL_INVALID_ARG_VALUE: return "invalid arg value";
    case CL_INVALID_ARG_SIZE: return "invalid arg size";
    case CL_INVALID_KERNEL_ARGS: return "invalid kernel args";
    case CL_INVALID_WORK_DIMENSION: return "invalid work dimension";
    case CL_INVALID_WORK_GROUP_SIZE: return "invalid work group size";
    case CL_INVALID_WORK_ITEM_SIZE: return "invalid work item size";
    case CL_INVALID_GLOBAL_OFFSET: return "invalid global offset";
    case CL_INVALID_EVENT_WAIT_LIST: return "invalid event wait list";
    case CL_INVALID_EVENT: return "invalid event";
    case CL_INVALID_OPERATION: return "invalid operation";
    case CL_INVALID_GL_OBJECT: return "invalid gl object";
    case CL_INVALID_BUFFER_SIZE: return "invalid buffer size";
    case CL_INVALID_MIP_LEVEL: return "invalid mip level";

#if defined(cl_khr_gl_sharing) && (cl_khr_gl_sharing >= 1)
    case CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR: return "invalid gl sharegroup reference number";
#endif

#ifdef CL_VERSION_1_1
    case CL_MISALIGNED_SUB_BUFFER_OFFSET: return "misaligned sub-buffer offset";
    case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST: return "exec status error for events in wait list";
    case CL_INVALID_GLOBAL_WORK_SIZE: return "invalid global work size";
#endif

    default: return "invalid/unknown error code";
  }
}

#define MAX_ENTRIES 10
#define BUF_SIZE 160000 // for log buffer
#define SOURCE_BUF_SIZE 16000

#define checkOpenclError(param) { openclAssert((param), __FILE__, __LINE__); }
void openclAssert(cl_int code, const char *file, int line)
{
        if (code != CL_SUCCESS) {
                fprintf(stderr, "openclAssert: error: %s, file %s, line %d\n", cl_error_to_str(code), file, line);
                exit(code);
        }
}

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

//extern __shared__ real shared_grid[];

void le_step_x(const int nx, const int ny, const real k1, const real k2, const le_material mat,
		const real* in_grid, real* grid)
{
	/*const int i = threadIdx.x + blockDim.x * blockIdx.x;
	const int j = threadIdx.y + blockDim.y * blockIdx.y;
	const int li = threadIdx.x;
	const int lj = threadIdx.y;

	real vx, vy, sxx, sxy, syy;
	le_w w_2, w_1, w, w1, w2, d;

	s_indx_vx(li, lj) = g_ind(in_grid, i, j, nx, ny, 0);
	s_indx_vy(li, lj) = g_ind(in_grid, i, j, nx, ny, 1);
	s_indx_sxx(li, lj) = g_ind(in_grid, i, j, nx, ny, 2);
	s_indx_sxy(li, lj) = g_ind(in_grid, i, j, nx, ny, 3);
	s_indx_syy(li, lj) = g_ind(in_grid, i, j, nx, ny, 4);
	if (li < 2) {
		s_indx_vx(li - 2, lj) = g_ind(in_grid, i - 2, j, nx, ny, 0);
		s_indx_vy(li - 2, lj) = g_ind(in_grid, i - 2, j, nx, ny, 1);
		s_indx_sxx(li - 2, lj) = g_ind(in_grid, i - 2, j, nx, ny, 2);
		s_indx_sxy(li - 2, lj) = g_ind(in_grid, i - 2, j, nx, ny, 3);
		s_indx_syy(li - 2, lj) = g_ind(in_grid, i - 2, j, nx, ny, 4);
	} else if (li >= blockDim.x - 2) {
		s_indx_vx(li + 2, lj) = g_ind(in_grid, i + 2, j, nx, ny, 0);
		s_indx_vy(li + 2, lj) = g_ind(in_grid, i + 2, j, nx, ny, 1);
		s_indx_sxx(li + 2, lj) = g_ind(in_grid, i + 2, j, nx, ny, 2);
		s_indx_sxy(li + 2, lj) = g_ind(in_grid, i + 2, j, nx, ny, 3);
		s_indx_syy(li + 2, lj) = g_ind(in_grid, i + 2, j, nx, ny, 4);
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
	ind_syy(i,j) = syy;*/
}

void le_step_y(const int nx, const int ny, const real k1, const real k2, const le_material mat,
		const real* in_grid, real* grid)
{
	/*const int i = threadIdx.x + blockDim.x * blockIdx.x;
	const int j = threadIdx.y + blockDim.y * blockIdx.y;
	const int li = threadIdx.x;
	const int lj = threadIdx.y;

	real vx, vy, sxx, sxy, syy;
	le_w w_2, w_1, w, w1, w2, d;

	s_indy_vx(li, lj) = g_ind(in_grid, i, j, nx, ny, 0);
	s_indy_vy(li, lj) = g_ind(in_grid, i, j, nx, ny, 1);
	s_indy_sxx(li, lj) = g_ind(in_grid, i, j, nx, ny, 2);
	s_indy_sxy(li, lj) = g_ind(in_grid, i, j, nx, ny, 3);
	s_indy_syy(li, lj) = g_ind(in_grid, i, j, nx, ny, 4);
	if (lj < 2) {
		s_indy_vx(li, lj - 2) = g_ind(in_grid, i, j - 2, nx, ny, 0);
		s_indy_vy(li, lj - 2) = g_ind(in_grid, i, j - 2, nx, ny, 1);
		s_indy_sxx(li, lj - 2) = g_ind(in_grid, i, j - 2, nx, ny, 2);
		s_indy_sxy(li, lj - 2) = g_ind(in_grid, i, j - 2, nx, ny, 3);
		s_indy_syy(li, lj - 2) = g_ind(in_grid, i, j - 2, nx, ny, 4);
	} else if (lj >= blockDim.y - 2) {
		s_indy_vx(li, lj + 2) = g_ind(in_grid, i, j + 2, nx, ny, 0);
		s_indy_vy(li, lj + 2) = g_ind(in_grid, i, j + 2, nx, ny, 1);
		s_indy_sxx(li, lj + 2) = g_ind(in_grid, i, j + 2, nx, ny, 2);
		s_indy_sxy(li, lj + 2) = g_ind(in_grid, i, j + 2, nx, ny, 3);
		s_indy_syy(li, lj + 2) = g_ind(in_grid, i, j + 2, nx, ny, 4);
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
	ind_syy(i,j) = syy;*/
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
	int i;
	
	const int nx = task->n.x;
	const int ny = task->n.y;
	
	const int grid_size = sizeof(real) * NODE_SIZE * nx * ny;

	const real k1x = task->dt * task->mat.c1 / task->h.x;
	const real k2x = task->dt * task->mat.c2 / task->h.x;
	const real k1y = task->dt * task->mat.c1 / task->h.y;
	const real k2y = task->dt * task->mat.c2 / task->h.y;

	// set sizes of blocks on gpu
	int block_dim_stepx_x = 256;
	int block_dim_stepx_y = 1;
	size_t local_work_size_x[2] = {block_dim_stepx_x, block_dim_stepx_y};
	size_t global_work_size_x[2] = {(nx + block_dim_stepx_x - 1) / block_dim_stepx_x * block_dim_stepx_x,
		(ny  + block_dim_stepx_y - 1) / block_dim_stepx_y * block_dim_stepx_y};
	//size_t global_work_size_x[2] = {(nx + block_dim_stepx_x - 1) / block_dim_stepx_x, (ny  + block_dim_stepx_y - 1) / block_dim_stepx_y};
	// dim3 threadsPerBlockX(block_dim_stepx_x, block_dim_stepx_y);
	// dim3 blocksPerGridX((nx + block_dim_stepx_x - 1) / block_dim_stepx_x, (ny  + block_dim_stepx_y - 1) / block_dim_stepx_y);
	int sharedMemSizeX = (block_dim_stepx_x + 4) * block_dim_stepx_y * NODE_SIZE * sizeof(real);
	
	int block_dim_stepy_x = 16;
	int block_dim_stepy_y = 16;
	size_t local_work_size_y[2] = {block_dim_stepy_x, block_dim_stepy_y};
	size_t global_work_size_y[2] = {(nx + block_dim_stepy_x - 1) / block_dim_stepy_x * block_dim_stepy_x, 
		(ny  + block_dim_stepy_y - 1) / block_dim_stepy_y * block_dim_stepy_y};
	//size_t global_work_size_y[2] = {(nx + block_dim_stepy_x - 1) / block_dim_stepy_x, (ny  + block_dim_stepy_y - 1) / block_dim_stepy_y};
	// dim3 threadsPerBlockY(block_dim_stepy_x, block_dim_stepy_y);
	// dim3 blocksPerGridY((nx + block_dim_stepy_x - 1) / block_dim_stepy_x, (ny  + block_dim_stepy_y - 1) / block_dim_stepy_y);
	int sharedMemSizeY = block_dim_stepy_x * (block_dim_stepy_y + 4) * NODE_SIZE * sizeof(real);
	
	double t;

	// variables to work with opencl
	char buf[BUF_SIZE];
	int err_code;
        cl_platform_id platforms[MAX_ENTRIES];
        cl_uint num_platforms;
	int platform_num = 0;
	cl_device_id devices[MAX_ENTRIES];
	cl_uint num_devices;
	int device_num = 0;
	cl_context_properties context_properties[3] = { CL_CONTEXT_PLATFORM, 0, 0} ;
	cl_context context;
	cl_command_queue command_queue;
	cl_program program;
	cl_build_status build_status;
	cl_kernel kernel_step_x, kernel_step_y;
	cl_mem d_grid1, d_grid2;
	char* source = malloc(sizeof(char) * SOURCE_BUF_SIZE);
	assert(source);
	
	// get platforms
        checkOpenclError(clGetPlatformIDs(MAX_ENTRIES, platforms, &num_platforms));
	if (num_platforms == 0) {
		printf("There are no OpenCL platforms, exiting...\n");
		exit(-1);
	}
	
        for (i = 0; i < num_platforms; i++)
        {
                printf("platform %d:\n", i);
                checkOpenclError(clGetPlatformInfo(platforms[i], CL_PLATFORM_PROFILE, BUF_SIZE * sizeof(char), buf, NULL));
                printf("\tbuf: %s\n", buf);
		checkOpenclError(clGetPlatformInfo(platforms[i], CL_PLATFORM_VERSION, BUF_SIZE * sizeof(char), buf, NULL));
                printf("\tversion: %s\n", buf);
		checkOpenclError(clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, BUF_SIZE * sizeof(char), buf, NULL));
                printf("\tname: %s\n", buf);
		checkOpenclError(clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR, BUF_SIZE * sizeof(char), buf, NULL));
                printf("\tvendor: %s\n", buf);
		checkOpenclError(clGetPlatformInfo(platforms[i], CL_PLATFORM_EXTENSIONS, BUF_SIZE * sizeof(char), buf, NULL));
                printf("\textensions: %s\n", buf);
        }
        
        if (num_platforms != 1) {
		printf("Choose platform number: ");
		assert(scanf("%d", &platform_num));
	}
        
        // get devices
	checkOpenclError(clGetDeviceIDs(platforms[platform_num], CL_DEVICE_TYPE_ALL, MAX_ENTRIES, devices, &num_devices));
	if (num_devices == 0) {
		printf("There are no OpenCL devices, exiting...\n");
		exit(-1);
	}
	
	for (i = 0; i < num_devices; i++)
	{
		checkOpenclError(clGetDeviceInfo(devices[i], CL_DEVICE_NAME, BUF_SIZE * sizeof(char), buf, NULL));
		printf("device %d: %s\n", i, buf);
	}
	
	if (num_devices != 1) {
		printf("Choose device number: ");
		assert(scanf("%d", &device_num));
	}
	
	// create context
	context_properties[1] = (cl_context_properties) platforms[platform_num];
	context = clCreateContext(context_properties, 1, &devices[device_num], NULL, NULL, &err_code);
	checkOpenclError(err_code);
	
	// create command queue
	command_queue = clCreateCommandQueue(context, devices[device_num], 0, &err_code);
	checkOpenclError(err_code);
	
	// allocate memory on gpu
	d_grid1 = clCreateBuffer(context, CL_MEM_READ_WRITE, grid_size, NULL, &err_code);
	checkOpenclError(err_code);
	d_grid2 = clCreateBuffer(context, CL_MEM_READ_WRITE, grid_size, NULL, &err_code);
	checkOpenclError(err_code);
	
	// copy grid to device
	checkOpenclError(clEnqueueWriteBuffer(command_queue, d_grid1, CL_TRUE, 0, grid_size, task->grid, 0, NULL, NULL));

	// load program source from file
	loadFromFile("kernel.cl", source, SOURCE_BUF_SIZE);
	
	// create program
	program = clCreateProgramWithSource(context, 1, (const char **) &source, NULL, &err_code);
	checkOpenclError(err_code);
	
	// build program
	err_code = clBuildProgram(program, 1, &devices[device_num], options, NULL, NULL);
	
	// print build info
	checkOpenclError(clGetProgramBuildInfo(program, devices[device_num], CL_PROGRAM_BUILD_LOG, BUF_SIZE * sizeof(char), buf, NULL));
	checkOpenclError(clGetProgramBuildInfo(program, devices[device_num], CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), &build_status, NULL));
	printf("build log:\n%s\n", buf);
	if (CL_BUILD_SUCCESS == build_status) printf("Build program: OK\n");
	else {
		printf("Build program: FAIL\n");
		exit(-1);
	}
	
	// check clBuildProgram
	checkOpenclError(err_code);
	
	// create kernel
	kernel_step_x = clCreateKernel(program, "le_step_x", &err_code);
	checkOpenclError(err_code);
	kernel_step_y = clCreateKernel(program, "le_step_y", &err_code);
	checkOpenclError(err_code);
	
	// set kernel args
	checkOpenclError(clSetKernelArg(kernel_step_x, 0, sizeof(nx), &nx));
	checkOpenclError(clSetKernelArg(kernel_step_x, 1, sizeof(ny), &ny));
	checkOpenclError(clSetKernelArg(kernel_step_x, 2, sizeof(k1x), &k1x));
	checkOpenclError(clSetKernelArg(kernel_step_x, 3, sizeof(k2x), &k2x));
	checkOpenclError(clSetKernelArg(kernel_step_x, 4, sizeof(task->mat), &(task->mat)));
	checkOpenclError(clSetKernelArg(kernel_step_x, 5, sizeof(d_grid1), &d_grid1));
	checkOpenclError(clSetKernelArg(kernel_step_x, 6, sizeof(d_grid2), &d_grid2));
	checkOpenclError(clSetKernelArg(kernel_step_x, 7, sharedMemSizeX, NULL));
	
	checkOpenclError(clSetKernelArg(kernel_step_y, 0, sizeof(nx), &nx));
	checkOpenclError(clSetKernelArg(kernel_step_y, 1, sizeof(ny), &ny));
	checkOpenclError(clSetKernelArg(kernel_step_y, 2, sizeof(k1x), &k1x));
	checkOpenclError(clSetKernelArg(kernel_step_y, 3, sizeof(k2x), &k2x));
	checkOpenclError(clSetKernelArg(kernel_step_y, 4, sizeof(task->mat), &(task->mat)));
	checkOpenclError(clSetKernelArg(kernel_step_y, 5, sizeof(d_grid2), &d_grid2));
	checkOpenclError(clSetKernelArg(kernel_step_y, 6, sizeof(d_grid1), &d_grid1));
	checkOpenclError(clSetKernelArg(kernel_step_y, 7, sharedMemSizeY, NULL));
	
	t = timer();
	// add kernel execution in queue
	for (i = 0; i < steps; i++)
	{
		checkOpenclError(clEnqueueNDRangeKernel(command_queue, kernel_step_x, 2, NULL, global_work_size_x, local_work_size_x, 0, NULL, NULL));
		checkOpenclError(clEnqueueNDRangeKernel(command_queue, kernel_step_y, 2, NULL, global_work_size_y, local_work_size_y, 0, NULL, NULL));
	}
	clFinish(command_queue);
	t = timer() - t;

	// return data to host
	checkOpenclError(clEnqueueReadBuffer(command_queue, d_grid1, CL_TRUE, 0, grid_size, task->grid, 0, NULL, NULL));

	// release kernel
	checkOpenclError(clReleaseKernel(kernel_step_x));
	checkOpenclError(clReleaseKernel(kernel_step_y));
	
	// release program
	checkOpenclError(clReleaseProgram(program));
	
	// release command queue
	checkOpenclError(clReleaseCommandQueue(command_queue));
	
	// release context
	checkOpenclError(clReleaseContext(context));
	
	free(source);
	
	return t;
}


