#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "sdecomp.h"
#include "memory.h"
#include "runge_kutta.h"
#include "domain.h"
#include "tdm.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "array_macros/domain/dxf.h"
#include "array_macros/domain/dxc.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/fluid/psi.h"

static const double g_pi = 3.14159265358979324;

// structure only used to solve Poisson equation
// NOTE: this is shared among normal and efficient solvers
//   thus some variables may not be used
typedef struct {
  bool is_initialised;
  void * restrict buf0;
  void * restrict buf1;
  fftw_plan fftw_plan_x[2];
  fftw_plan fftw_plan_y[2];
  size_t tdm_sizes[2];
  tdm_info_t * tdm_info;
  double * evals;
  sdecomp_transpose_plan_t * r_transposer_x1_to_y1;
  sdecomp_transpose_plan_t * r_transposer_y1_to_x1;
  sdecomp_transpose_plan_t * c_transposer_x1_to_y1;
  sdecomp_transpose_plan_t * c_transposer_y1_to_x1;
} poisson_solver_t;

/* initialise Poisson solver */
// several pencils for different data types are treated
//   and thus this source is very complicated
// used prefixes are as follows:
//   r_: real    (double)       type
//   c_: complex (fftw_complex) type
//
//   gl_    : global array size (not pencils)
//   x1pncl_: each pencl (x1, y1, ...)

/* size of domain and pencils */
// NOTE: define globally to reduce the number of arguments of static functions
// global domain size in real space
static size_t r_gl_sizes[NDIMS] = {0};
// global domain size in complex space
static size_t c_gl_sizes[NDIMS] = {0};
// local domain size (x1 pencil) in real space
static size_t r_x1pncl_sizes[NDIMS] = {0};
// local domain size (y1 pencil) in real space
static size_t r_y1pncl_sizes[NDIMS] = {0};
// local domain size (y1 pencil) in complex space
static size_t c_y1pncl_sizes[NDIMS] = {0};
// local domain size (x1 pencil) in complex space
static size_t c_x1pncl_sizes[NDIMS] = {0};

static size_t prod(
    const size_t sizes[NDIMS]
){
  // compute the product of the given vector
  size_t nitems = 1;
  for(size_t dim = 0; dim < NDIMS; dim++){
    nitems *= sizes[dim];
  }
  return nitems;
}

static int report_failure(
    const char type[]
){
  // function to just dump error message and abort
  FILE * stream = stderr;
  fprintf(stream, "Poisson solver, initialisation failed: %s\n", type);
  fprintf(stream, "  FFTW:    A possible reason is you link Intel-MKL lib\n");
  fprintf(stream, "           Make sure you use FFTW3 directly,\n");
  fprintf(stream, "           NOT its crazy wrapper offered by MKL\n");
  fprintf(stream, "  SDECOMP: Check sdecomp.log and check arguments\n");
  fprintf(stream, "           If they are all correct, PLEASE CONTACT ME\n");
  fflush(stream);
  return 0;
}

static size_t max(
    const size_t val0,
    const size_t val1
){
  if(val0 > val1){
    return val0;
  }else{
    return val1;
  }
}

static int compute_pencil_sizes(
    const domain_t * domain
){
  // NOTE: those variables are defined globally at the top of this file
  //   to reduce the nhumber of arguments which functions take
  // global domain size in real space
  const sdecomp_info_t * info = domain->info;
  r_gl_sizes[0] = domain->glsizes[0];
  r_gl_sizes[1] = domain->glsizes[1];
  // global domain size in complex space
  // NOTE: Hermite symmetry in y
  c_gl_sizes[0] = domain->glsizes[0];
  c_gl_sizes[1] = domain->glsizes[1] / 2 + 1;
  // local domain sizes
  for(sdecomp_dir_t dim = 0; dim < NDIMS; dim++){
    if(0 != sdecomp.get_pencil_mysize(info, SDECOMP_X1PENCIL, dim, r_gl_sizes[dim], r_x1pncl_sizes + dim)) return 1;
    if(0 != sdecomp.get_pencil_mysize(info, SDECOMP_Y1PENCIL, dim, r_gl_sizes[dim], r_y1pncl_sizes + dim)) return 1;
    if(0 != sdecomp.get_pencil_mysize(info, SDECOMP_Y1PENCIL, dim, c_gl_sizes[dim], c_y1pncl_sizes + dim)) return 1;
    if(0 != sdecomp.get_pencil_mysize(info, SDECOMP_X1PENCIL, dim, c_gl_sizes[dim], c_x1pncl_sizes + dim)) return 1;
  }
  return 0;
}

static int allocate_buffers(
    poisson_solver_t * poisson_solver
){
  // although there are bunch of pencils involved,
  //   two buffers are enough to do the job,
  //   which are allocated here
  void * restrict * buf0 = &poisson_solver->buf0;
  void * restrict * buf1 = &poisson_solver->buf1;
  const size_t r_dsize = sizeof(double);
  const size_t c_dsize = sizeof(fftw_complex);
  size_t buf0_bytes = 0;
  size_t buf1_bytes = 0;
  // r_x1pncl -> rotate -> r_y1pncl -> FFT -> c_y1pncl -> rotate -> c_x1pncl
  // buffer0               buffer1            buffer0               buffer1
  buf0_bytes = max(buf0_bytes, r_dsize * prod(r_x1pncl_sizes));
  buf0_bytes = max(buf0_bytes, c_dsize * prod(c_y1pncl_sizes));
  buf1_bytes = max(buf1_bytes, r_dsize * prod(r_y1pncl_sizes));
  buf1_bytes = max(buf1_bytes, c_dsize * prod(c_x1pncl_sizes));
  // allocate them using fftw_malloc to enforce them 16bit-aligned for SIMD
  *buf0 = fftw_malloc(buf0_bytes);
  if(NULL == *buf0){
    fprintf(stderr, "FATAL: fftw_malloc failed (requested %zu bytes)\n", buf0_bytes);
    fflush(stderr);
    return 1;
  }
  *buf1 = fftw_malloc(buf1_bytes);
  if(NULL == *buf1){
    fprintf(stderr, "FATAL: fftw_malloc failed (requested %zu bytes)\n", buf1_bytes);
    fflush(stderr);
    return 1;
  }
  return 0;
}

static int init_tri_diagonal_solver(
    const domain_t * domain,
    poisson_solver_t * poisson_solver
){
  // N x N tri-diagonal matrix,
  //   which are solved for M times
  // tdm_sizes[0] = N, tdm_sizes[1] = M
  // since lower- and upper-diagonal components are
  //   independent to y, z directions and time,
  //   we compute here and re-use them
  // center-diagonal components are, on the other hand,
  //   dependent on time and thus needs to compute everytime
  //   in the solver
  size_t * restrict tdm_sizes = poisson_solver->tdm_sizes;
  tdm_info_t ** tdm_info = &poisson_solver->tdm_info;
  // in x: d^2p / dx^2 = q
  tdm_sizes[0] = c_x1pncl_sizes[0];
  tdm_sizes[1] = c_x1pncl_sizes[1];
  if(0 != tdm.construct(
    /* size of system */ tdm_sizes[0],
    /* number of rhs  */ 1,
    /* is periodic    */ false,
    /* is complex     */ true,
    /* output         */ tdm_info
  )) return 1;
  // initialise tri-diagonal matrix in x direction
  double * tdm_l = NULL;
  double * tdm_u = NULL;
  tdm.get_l(*tdm_info, &tdm_l);
  tdm.get_u(*tdm_info, &tdm_u);
  const double * dxf = domain->dxf;
  const double * dxc = domain->dxc;
  for(size_t i = 1; i <= tdm_sizes[0]; i++){
    // N.B. loop from i = 1 to use DXC and DXF macros,
    //   which should be assigned to tdm_[lu] at i = 0
    tdm_l[i-1] = 1. / DXC(i  ) / DXF(i  );
    tdm_u[i-1] = 1. / DXC(i+1) / DXF(i  );
  }
  return 0;
}

static int init_pencil_rotations(
    const domain_t * domain,
    poisson_solver_t * poisson_solver
){
  const sdecomp_info_t * info = domain->info;
  const size_t r_dsize = sizeof(double);
  const size_t c_dsize = sizeof(fftw_complex);
  if(0 != sdecomp.transpose.construct(info, SDECOMP_X1PENCIL, SDECOMP_Y1PENCIL, r_gl_sizes, r_dsize, &poisson_solver->r_transposer_x1_to_y1)){
    report_failure("SDECOMP x1 to y1 for real");
    return 1;
  }
  if(0 != sdecomp.transpose.construct(info, SDECOMP_Y1PENCIL, SDECOMP_X1PENCIL, r_gl_sizes, r_dsize, &poisson_solver->r_transposer_y1_to_x1)){
    report_failure("SDECOMP y1 to x1 for real");
    return 1;
  }
  if(0 != sdecomp.transpose.construct(info, SDECOMP_X1PENCIL, SDECOMP_Y1PENCIL, c_gl_sizes, c_dsize, &poisson_solver->c_transposer_x1_to_y1)){
    report_failure("SDECOMP x1 to y1 for complex");
    return 1;
  }
  if(0 != sdecomp.transpose.construct(info, SDECOMP_Y1PENCIL, SDECOMP_X1PENCIL, c_gl_sizes, c_dsize, &poisson_solver->c_transposer_y1_to_x1)){
    report_failure("SDECOMP y1 to x1 for complex");
    return 1;
  }
  return 0;
}

static int init_ffts(
    poisson_solver_t * poisson_solver
){
  const unsigned flags = FFTW_PATIENT | FFTW_DESTROY_INPUT;
  // NOTE: two buffers should be properly given
  //   see "allocate_buffers" above
  // y, real / complex
  {
    fftw_plan * fplan = &poisson_solver->fftw_plan_y[0];
    fftw_plan * bplan = &poisson_solver->fftw_plan_y[1];
    const int r_signal_length = r_y1pncl_sizes[SDECOMP_YDIR];
    const int c_signal_length = c_y1pncl_sizes[SDECOMP_YDIR];
    const int repeat_for = r_y1pncl_sizes[SDECOMP_XDIR];
    *fplan = fftw_plan_many_dft_r2c(
        1, &r_signal_length, repeat_for,
        poisson_solver->buf1, NULL, 1, r_signal_length,
        poisson_solver->buf0, NULL, 1, c_signal_length,
        flags
    );
    *bplan = fftw_plan_many_dft_c2r(
        1, &r_signal_length, repeat_for,
        poisson_solver->buf0, NULL, 1, c_signal_length,
        poisson_solver->buf1, NULL, 1, r_signal_length,
        flags
    );
    if(NULL == *fplan){
      report_failure("FFTW y-forward");
      return 1;
    }
    if(NULL == *bplan){
      report_failure("FFTW y-backward");
      return 1;
    }
  }
  return 0;
}

static int init_eigenvalues(
    const domain_t * domain,
    poisson_solver_t * poisson_solver
){
  const sdecomp_info_t * info = domain->info;
  double ** evals = &poisson_solver->evals;
  // x1 pencil, DFT in y
  const sdecomp_pencil_t pencil = SDECOMP_X1PENCIL;
  const double signal_lengths[NDIMS - 1] = {
    r_gl_sizes[SDECOMP_YDIR],
  };
  size_t mysizes[NDIMS - 1] = {0};
  sdecomp.get_pencil_mysize(info, pencil, SDECOMP_YDIR, c_gl_sizes[SDECOMP_YDIR], mysizes);
  size_t offsets[NDIMS - 1] = {0};
  sdecomp.get_pencil_offset(info, pencil, SDECOMP_YDIR, c_gl_sizes[SDECOMP_YDIR], offsets);
  const double gridsizes[NDIMS - 1] = {
    domain->lengths[SDECOMP_YDIR] / r_gl_sizes[SDECOMP_YDIR],
  };
  // initialise eigenvalues in homogeneous directions
  *evals = memory_calloc(mysizes[0], sizeof(double));
  for(size_t cnt = 0, j = offsets[0]; j < mysizes[0] + offsets[0]; j++, cnt++){
    (*evals)[cnt] =
      - 4. / pow(gridsizes[0], 2.) * pow(
        sin( g_pi * j / signal_lengths[0] ),
        2.
    );
  }
  return 0;
}

static int init_poisson_solver(
    const domain_t * domain,
    poisson_solver_t * poisson_solver
){
  // check domain size (global, local, pencils)
  if(0 != compute_pencil_sizes(domain)) return 1;
  // initialise each part of poisson_solver_t
  if(0 != allocate_buffers(poisson_solver))                 return 1;
  if(0 != init_tri_diagonal_solver(domain, poisson_solver)) return 1;
  if(0 != init_pencil_rotations(domain, poisson_solver))    return 1;
  if(0 != init_ffts(poisson_solver))                        return 1;
  if(0 != init_eigenvalues(domain, poisson_solver))         return 1;
  poisson_solver->is_initialised = true;
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if(root == myrank){
    printf("DFT-based solver is used\n");
  }
  return 0;
}

static int assign_input(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    const fluid_t * fluid,
    double * restrict rhs
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict dxf = domain->dxf;
  const double dy = domain->dy;
  const double * restrict ux = fluid->ux.data;
  const double * restrict uy = fluid->uy.data;
  // normalise FFT beforehand
  const double norm = 1. * domain->glsizes[1];
  const double prefactor = 1. / (rkcoefs[rkstep][rk_g] * dt) / norm;
  for(int cnt = 0, j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++, cnt++){
      const double dx = DXF(i  );
      const double ux_xm = UX(i  , j  );
      const double ux_xp = UX(i+1, j  );
      const double uy_ym = UY(i  , j  );
      const double uy_yp = UY(i  , j+1);
      rhs[cnt] = prefactor * (
         + (ux_xp - ux_xm) / dx
         + (uy_yp - uy_ym) / dy
      );
    }
  }
  return 0;
}

static int extract_output(
    const domain_t * domain,
    const double * restrict rhs,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  double * restrict psi = fluid->psi.data;
  for(int cnt = 0, j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++, cnt++){
      PSI(i, j) = rhs[cnt];
    }
  }
  fluid_update_boundaries_psi(domain, &fluid->psi);
  return 0;
}

static int solve_linear_systems(
    poisson_solver_t * poisson_solver
){
  // size of system (length) and how many such systems to be solved
  // NOTE: although size_of_system is the same as tdm.get_size gives,
  //   repeat_for is different from what tdm.get_nrhs returns (=1)
  //   here repeat_for is the degree of freedom in the wavespace
  const size_t size_of_system = poisson_solver->tdm_sizes[0];
  const size_t repeat_for     = poisson_solver->tdm_sizes[1];
  // tri-diagonal matrix
  tdm_info_t * tdm_info = poisson_solver->tdm_info;
  double * restrict tdm_l = NULL;
  double * restrict tdm_u = NULL;
  double * restrict tdm_c = NULL;
  tdm.get_l(tdm_info, &tdm_l);
  tdm.get_u(tdm_info, &tdm_u);
  tdm.get_c(tdm_info, &tdm_c);
  // eigenvalues coming from Fourier projection
  const double * restrict evals = poisson_solver->evals;
  fftw_complex * restrict rhs = poisson_solver->buf1;
  for(size_t m = 0; m < repeat_for; m++){
    // set center diagonal components
    for(size_t n = 0; n < size_of_system; n++){
      tdm_c[n] = - tdm_l[n] - tdm_u[n] + evals[m];
    }
    // boundary treatment (Neumann boundary condition)
    tdm_c[               0] += tdm_l[               0];
    tdm_c[size_of_system-1] += tdm_u[size_of_system-1];
    tdm.solve(tdm_info, rhs + m * size_of_system);
  }
  return 0;
}

/**
 * @brief compute scalar potential psi to correct velocity
 * @param[in]     domain : information about domain decomposition and size
 * @param[in]     rkstep : Runge-Kutta step
 * @param[in]     dt     : time step size
 * @param[in,out] fluid  : velocity (in), scalar potential psi (out)
 * @return               : (success) 0
 *                       : (failure) 1
 */
int fluid_compute_potential_dft(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
){
  static poisson_solver_t poisson_solver = {
    .is_initialised = false,
  };
  // initialise Poisson solver
  if(!poisson_solver.is_initialised){
    if(0 != init_poisson_solver(domain, &poisson_solver)){
      // failed to initialise Poisson solver
      return 1;
    }
  }
  // compute right-hand side of Poisson equation
  // assigned to buf0
  assign_input(domain, rkstep, dt, fluid, poisson_solver.buf0);
  // solve the equation
  // transpose real x1pencil to y1pencil
  // from buf0 to buf1
  sdecomp.transpose.execute(
      poisson_solver.r_transposer_x1_to_y1,
      poisson_solver.buf0,
      poisson_solver.buf1
  );
  // project y to wave space
  // f(x, y)    -> f(x, k_y)
  // f(x, y, z) -> f(x, k_y, z)
  // from buf1 to buf0
  fftw_execute(poisson_solver.fftw_plan_y[0]);
  // transpose complex y1pencil to x1pencil
  // from buf0 to buf1
  sdecomp.transpose.execute(
      poisson_solver.c_transposer_y1_to_x1,
      poisson_solver.buf0,
      poisson_solver.buf1
  );
  // solve linear systems
  solve_linear_systems(&poisson_solver);
  // transpose complex x1pencil to y1pencil
  // from buf1 to buf0
  sdecomp.transpose.execute(
      poisson_solver.c_transposer_x1_to_y1,
      poisson_solver.buf1,
      poisson_solver.buf0
  );
  // project y to physical space
  // f(x, k_y)    -> f(x, y)
  // f(x, k_y, z) -> f(x, y, z)
  // from buf0 to buf1
  fftw_execute(poisson_solver.fftw_plan_y[1]);
  // transpose real y1pencil to x1pencil
  // from buf1 to buf0
  sdecomp.transpose.execute(
      poisson_solver.r_transposer_y1_to_x1,
      poisson_solver.buf1,
      poisson_solver.buf0
  );
  extract_output(domain, poisson_solver.buf0, fluid);
  return 0;
}

