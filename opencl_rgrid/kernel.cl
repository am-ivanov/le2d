#include "le_types.h"

#define s_indx_vx(i, j) (*(shared_grid + ((i)+2) + (j) * (get_local_size(0) + 4) + 0 * (get_local_size(0) + 4) * get_local_size(1)))
#define s_indx_vy(i, j) (*(shared_grid + ((i)+2) + (j) * (get_local_size(0) + 4) + 1 * (get_local_size(0) + 4) * get_local_size(1)))
#define s_indx_sxx(i, j) (*(shared_grid + ((i)+2) + (j) * (get_local_size(0) + 4) + 2 * (get_local_size(0) + 4) * get_local_size(1)))
#define s_indx_sxy(i, j) (*(shared_grid + ((i)+2) + (j) * (get_local_size(0) + 4) + 3 * (get_local_size(0) + 4) * get_local_size(1)))
#define s_indx_syy(i, j) (*(shared_grid + ((i)+2) + (j) * (get_local_size(0) + 4) + 4 * (get_local_size(0) + 4) * get_local_size(1)))

#define s_indy_vx(i, j) (*(shared_grid + (i) + ((j)+2) * get_local_size(0) + 0 * get_local_size(0) * (get_local_size(1) + 4)))
#define s_indy_vy(i, j) (*(shared_grid + (i) + ((j)+2) * get_local_size(0) + 1 * get_local_size(0) * (get_local_size(1) + 4)))
#define s_indy_sxx(i, j) (*(shared_grid + (i) + ((j)+2) * get_local_size(0) + 2 * get_local_size(0) * (get_local_size(1) + 4)))
#define s_indy_sxy(i, j) (*(shared_grid + (i) + ((j)+2) * get_local_size(0) + 3 * get_local_size(0) * (get_local_size(1) + 4)))
#define s_indy_syy(i, j) (*(shared_grid + (i) + ((j)+2) * get_local_size(0) + 4 * get_local_size(0) * (get_local_size(1) + 4)))

#define vnorm(v) (sqrt(v.x * v.x + v.y * v.y))

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

inline real le_min(real a, real b) { return a > b ? b : a; }
inline real le_max(real a, real b) { return a > b ? a : b; }
inline real le_max3(real a, real b, real c) { return le_max(a, le_max(b, c)); }

#define limiter_minmod(r) (le_max(CONST_ZERO, le_min(CONST_ONE, (r))))
#define limiter_cir(r) (CONST_ZERO)
#define limiter_superbee(r) (le_max3(CONST_ZERO, le_min(CONST_ONE, CONST_TWO * r), le_min(CONST_TWO, r)))

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

	const real k = CONST_HALF * (CONST_ONE - c);

	return c * ((f_12 - f12) * k - r1);
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

inline void inc_y(const le_material *m, const real vx, const real vy, const real sxx, const real sxy, const real syy,
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

inline void reconstruct(const le_w ppu, const le_w pu, const le_w u, const le_w nu, const le_w nnu, const real k1, const real k2, le_w *d)
{
	d->w1 = tvd2(k1, ppu.w1, pu.w1, u.w1, nu.w1); // c1
	d->w2 = tvd2(k1, nnu.w2, nu.w2, u.w2, pu.w2); // -c1
	d->w3 = tvd2(k2, ppu.w3, pu.w3, u.w3, nu.w3); // c2
	d->w4 = tvd2(k2, nnu.w4, nu.w4, u.w4, pu.w4); // -c2
}

inline real g_ind_x(__global const real* grid, int i, int j, const int nx, const int ny, const int node)
{
	if (i < 0) i = 0;
	if (i >= nx) i = nx - 1;
	return (*(grid + (i+2) + (j+2) * (nx+4) + node * (nx+4) * (ny+4)));
}

inline real g_ind_y(__global const real* grid, int i, int j, const int nx, const int ny, const int node, const unsigned side_edge)
{
	if (j < 0 && (side_edge & 1u)) j = 0;
	if (j >= ny && (side_edge & 2u)) j = ny - 1;
	return (*(grid + (i+2) + (j+2) * (nx+4) + node * (nx+4) * (ny+4)));
}

#define g_ind(g, i, j, nx, ny, node) \
	(g)[(i + 2) + (j + 2) * (nx + 4) + (node) * (nx + 4) * (ny + 4)]

__kernel void le_step_x(const int nx, const int ny, const real k1, const real k2, const le_material mat,
		__global const real* in_grid, __global real* grid, __local real* shared_grid)
{
	const int i = get_global_id(0);
	const int j = get_global_id(1);
	const int li = get_local_id(0);
	const int lj = get_local_id(1);

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
	} else if (li >= get_local_size(0) - 2) {
		s_indx_vx(li + 2, lj) = g_ind_x(in_grid, i + 2, j, nx, ny, 0);
		s_indx_vy(li + 2, lj) = g_ind_x(in_grid, i + 2, j, nx, ny, 1);
		s_indx_sxx(li + 2, lj) = g_ind_x(in_grid, i + 2, j, nx, ny, 2);
		s_indx_sxy(li + 2, lj) = g_ind_x(in_grid, i + 2, j, nx, ny, 3);
		s_indx_syy(li + 2, lj) = g_ind_x(in_grid, i + 2, j, nx, ny, 4);
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	if (i >= nx || j >= ny) return;

	omega_x(&mat, s_indx_vx(li-2, lj), s_indx_vy(li-2, lj), s_indx_sxx(li-2, lj), s_indx_sxy(li-2, lj), &w_2);
	omega_x(&mat, s_indx_vx(li-1, lj), s_indx_vy(li-1, lj), s_indx_sxx(li-1, lj), s_indx_sxy(li-1, lj), &w_1);
	omega_x(&mat, s_indx_vx(li, lj), s_indx_vy(li, lj), s_indx_sxx(li, lj), s_indx_sxy(li, lj), &w);
	omega_x(&mat, s_indx_vx(li+1, lj), s_indx_vy(li+1, lj), s_indx_sxx(li+1, lj), s_indx_sxy(li+1, lj), &w1);
	omega_x(&mat, s_indx_vx(li+2, lj), s_indx_vy(li+2, lj), s_indx_sxx(li+2, lj), s_indx_sxy(li+2, lj), &w2);

	reconstruct(w_2, w_1, w, w1, w2, k1, k2, &d);
	inc_x(&mat, s_indx_vx(li,lj), s_indx_vy(li,lj), s_indx_sxx(li,lj), s_indx_sxy(li,lj), s_indx_syy(li,lj), &vx, &vy, &sxx, &sxy, &syy, &d);

	g_ind(grid, i, j, nx, ny, 0) = vx;
	g_ind(grid, i, j, nx, ny, 1) = vy;
	g_ind(grid, i, j, nx, ny, 2) = sxx;
	g_ind(grid, i, j, nx, ny, 3) = sxy;
	g_ind(grid, i, j, nx, ny, 4) = syy;
}

__kernel void le_step_y(const int nx, const int ny, const real k1, const real k2, const le_material mat,
		__global const real* in_grid, __global real* grid, __local real* shared_grid, const unsigned side_edge)
{
	const int i = get_global_id(0);
	const int j = get_global_id(1);
	const int li = get_local_id(0);
	const int lj = get_local_id(1);

	real vx, vy, sxx, sxy, syy;
	le_w w_2, w_1, w, w1, w2, d;

	s_indy_vx(li, lj) = g_ind_y(in_grid, i, j, nx, ny, 0, side_edge);
	s_indy_vy(li, lj) = g_ind_y(in_grid, i, j, nx, ny, 1, side_edge);
	s_indy_sxx(li, lj) = g_ind_y(in_grid, i, j, nx, ny, 2, side_edge);
	s_indy_sxy(li, lj) = g_ind_y(in_grid, i, j, nx, ny, 3, side_edge);
	s_indy_syy(li, lj) = g_ind_y(in_grid, i, j, nx, ny, 4, side_edge);
	if (lj < 2) {
		s_indy_vx(li, lj - 2) = g_ind_y(in_grid, i, j - 2, nx, ny, 0, side_edge);
		s_indy_vy(li, lj - 2) = g_ind_y(in_grid, i, j - 2, nx, ny, 1, side_edge);
		s_indy_sxx(li, lj - 2) = g_ind_y(in_grid, i, j - 2, nx, ny, 2, side_edge);
		s_indy_sxy(li, lj - 2) = g_ind_y(in_grid, i, j - 2, nx, ny, 3, side_edge);
		s_indy_syy(li, lj - 2) = g_ind_y(in_grid, i, j - 2, nx, ny, 4, side_edge);
	} else if (lj >= get_local_size(1) - 2) {
		s_indy_vx(li, lj + 2) = g_ind_y(in_grid, i, j + 2, nx, ny, 0, side_edge);
		s_indy_vy(li, lj + 2) = g_ind_y(in_grid, i, j + 2, nx, ny, 1, side_edge);
		s_indy_sxx(li, lj + 2) = g_ind_y(in_grid, i, j + 2, nx, ny, 2, side_edge);
		s_indy_sxy(li, lj + 2) = g_ind_y(in_grid, i, j + 2, nx, ny, 3, side_edge);
		s_indy_syy(li, lj + 2) = g_ind_y(in_grid, i, j + 2, nx, ny, 4, side_edge);
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	if (i >= nx || j >= ny) return;

	omega_y(&mat, s_indy_vx(li, lj-2), s_indy_vy(li, lj-2), s_indy_sxy(li, lj-2), s_indy_syy(li, lj-2), &w_2);
	omega_y(&mat, s_indy_vx(li, lj-1), s_indy_vy(li, lj-1), s_indy_sxy(li, lj-1), s_indy_syy(li, lj-1), &w_1);
	omega_y(&mat, s_indy_vx(li, lj), s_indy_vy(li, lj), s_indy_sxy(li, lj), s_indy_syy(li, lj), &w);
	omega_y(&mat, s_indy_vx(li, lj+1), s_indy_vy(li, lj+1), s_indy_sxy(li, lj+1), s_indy_syy(li, lj+1), &w1);
	omega_y(&mat, s_indy_vx(li, lj+2), s_indy_vy(li, lj+2), s_indy_sxy(li, lj+2), s_indy_syy(li, lj+2), &w2);

	reconstruct(w_2, w_1, w, w1, w2, k1, k2, &d);
	inc_y(&mat, s_indy_vx(li,lj), s_indy_vy(li,lj), s_indy_sxx(li,lj), s_indy_sxy(li,lj), s_indy_syy(li,lj), &vx, &vy, &sxx, &sxy, &syy, &d);

	g_ind(grid, i, j, nx, ny, 0) = vx;
	g_ind(grid, i, j, nx, ny, 1) = vy;
	g_ind(grid, i, j, nx, ny, 2) = sxx;
	g_ind(grid, i, j, nx, ny, 3) = sxy;
	g_ind(grid, i, j, nx, ny, 4) = syy;
}