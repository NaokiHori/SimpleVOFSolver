#include <math.h>
#include <string.h>
#include "runge_kutta.h"
#include "domain.h"
#include "fluid.h"
#include "interface.h"
#include "interface_solver.h"
#include "../internal.h"
#include "internal.h"
#include "array_macros/domain/dxf.h"
#include "array_macros/interface/vof.h"
#include "array_macros/interface/flxx.h"
#include "array_macros/interface/flxy.h"
#if NDIMS == 3
#include "array_macros/interface/flxz.h"
#endif
#include "array_macros/interface/src.h"

double indicator(
    const normal_t n,
    const vector_t x
){
  // diffused surface representation
  return 1. / (1. + exp(-2. * vofbeta * (
          + n[0] * x[0]
          + n[1] * x[1]
#if NDIMS == 3
          + n[2] * x[2]
#endif
          + n[NDIMS]
  )));
}

static int reset_srcs(
    const size_t rkstep,
    array_t * restrict srca,
    array_t * restrict srcb
){
  // copy previous k-step source term and reset
  if(0 != rkstep){
    // stash previous RK source term,
    //   which is achieved by swapping
    //   the pointers to "data"
    double * tmp = srca->data;
    srca->data = srcb->data;
    srcb->data = tmp;
  }
  // zero-clear current RK source terms
  memset(srca->data, 0, srca->datasize);
  return 0;
}

static int interface_compute_rhs(
    const domain_t * domain,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict dxf = domain->dxf;
  const double            dy  = domain->dy;
#if NDIMS == 3
  const double            dz  = domain->dz;
#endif
  const double * restrict flxx = interface->flxx.data;
  const double * restrict flxy = interface->flxy.data;
#if NDIMS == 3
  const double * restrict flxz = interface->flxz.data;
#endif
  double * restrict src = interface->src[rk_a].data;
  // compute right-hand-side of advection equation, x flux | 23
#if NDIMS == 2
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      const double dx = DXF(i  );
      SRC(i, j) += 1. / dx * (
          + FLXX(i  , j  )
          - FLXX(i+1, j  )
      );
    }
  }
#else
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        const double dx = DXF(i  );
        SRC(i, j, k) += 1. / dx * (
            + FLXX(i  , j  , k  )
            - FLXX(i+1, j  , k  )
        );
      }
    }
  }
#endif
  // compute right-hand-side of advection equation, y flux | 21
#if NDIMS == 2
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      SRC(i, j) += 1. / dy * (
          + FLXY(i  , j  )
          - FLXY(i  , j+1)
      );
    }
  }
#else
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        SRC(i, j, k) += 1. / dy * (
            + FLXY(i  , j  , k  )
            - FLXY(i  , j+1, k  )
        );
      }
    }
  }
#endif
#if NDIMS == 3
  // compute right-hand-side of advection equation, z flux | 10
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        SRC(i, j, k) += 1. / dz * (
            + FLXZ(i  , j  , k  )
            - FLXZ(i  , j  , k+1)
        );
      }
    }
  }
#endif
  return 0;
}

static int interface_advect_vof(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  double * restrict vof = interface->vof.data;
  // update vof, alpha contribution | 19
  {
    const double coef = rkcoefs[rkstep][rk_a];
    const double * restrict src = interface->src[rk_a].data;
#if NDIMS == 2
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        VOF(i, j) += dt * coef * SRC(i, j);
      }
    }
#else
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 1; i <= isize; i++){
          VOF(i, j, k) += dt * coef * SRC(i, j, k);
        }
      }
    }
#endif
  }
  // update vof, beta contribution | 19
  if(0 != rkstep){
    const double coef = rkcoefs[rkstep][rk_b];
    const double * restrict src = interface->src[rk_b].data;
#if NDIMS == 2
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        VOF(i, j) += dt * coef * SRC(i, j);
      }
    }
#else
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 1; i <= isize; i++){
          VOF(i, j, k) += dt * coef * SRC(i, j, k);
        }
      }
    }
#endif
  }
  return 0;
}

int interface_update_vof(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    const fluid_t * fluid,
    interface_t * interface
){
  reset_srcs(rkstep, interface->src + rk_a, interface->src + rk_b);
  compute_flux_x(domain, fluid, interface);
  compute_flux_y(domain, fluid, interface);
#if NDIMS == 3
  compute_flux_z(domain, fluid, interface);
#endif
  interface_compute_rhs(domain, interface);
  interface_advect_vof(domain, rkstep, dt, interface);
  interface_update_boundaries_vof(domain, &interface->vof);
  return 0;
}

