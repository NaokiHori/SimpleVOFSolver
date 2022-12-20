#include <string.h>
#include "common.h"
#include "domain.h"
#include "fluid.h"
#include "interface.h"
#include "arrays/fluid.h"
#include "arrays/interface.h"
#include "internal.h"



static int interface_compute_flux_x(const domain_t *domain, const fluid_t *fluid, interface_t *interface){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict ux = fluid->ux;
  const double * restrict gps = interface->gps;
  const double * restrict gws = interface->gws;
  const double * restrict vof = interface->vof;
  const nrml_t * restrict normal = interface->normal;
  double * restrict voffluxx = interface->voffluxx;
  for(int j = 1; j <= jsize; j++){
    for(int i = 2; i <= isize; i++){
      /* ! use upwind information ! 9 ! */
      int ii;
      double x;
      if(UX(i, j) < 0.){
        ii = i;
        x = -0.5;
      }else{
        ii = i-1;
        x = +0.5;
      }
      /* ! evaluate flux ! 12 ! */
      double flux = 0.;
      if(VOF(ii, j) < VOFMIN || 1.-VOFMIN < VOF(ii, j)){
        flux = VOF(ii, j);
      }else{
        double a10 = NORMAL(ii, j).a10;
        double a01 = NORMAL(ii, j).a01;
        double a00 = NORMAL(ii, j).a00;
        for(int jj = 0; jj < ORDER_GAUSS; jj++){
          flux += gws[jj] * H(a00, a10, a01, x, gps[jj]);
        }
      }
      VOFFLUXX(i, j) = flux * UX(i, j);
    }
  }
  return 0;
}

static int interface_compute_flux_y(const domain_t *domain, const fluid_t *fluid, interface_t *interface){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict uy = fluid->uy;
  const double * restrict gps = interface->gps;
  const double * restrict gws = interface->gws;
  const double * restrict vof = interface->vof;
  const nrml_t * restrict normal = interface->normal;
  double * restrict voffluxy = interface->voffluxy;
  for(int j = 1; j <= jsize+1; j++){
    for(int i = 1; i <= isize; i++){
      /* ! use upwind information ! 9 ! */
      int jj;
      double y;
      if(UY(i, j) < 0.){
        jj = j;
        y = -0.5;
      }else{
        jj = j-1;
        y = +0.5;
      }
      /* ! evaluate flux ! 12 ! */
      double flux = 0.;
      if(VOF(i, jj) < VOFMIN || 1.-VOFMIN < VOF(i, jj)){
        flux = VOF(i, jj);
      }else{
        double a10 = NORMAL(i, jj).a10;
        double a01 = NORMAL(i, jj).a01;
        double a00 = NORMAL(i, jj).a00;
        for(int ii = 0; ii < ORDER_GAUSS; ii++){
          flux += gws[ii] * H(a00, a10, a01, gps[ii], y);
        }
      }
      VOFFLUXY(i, j) = flux * UY(i, j);
    }
  }
  return 0;
}

static int interface_compute_rhs(const domain_t *domain, interface_t *interface){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict dxf = domain->dxf;
  const double            dy  = domain->dy;
  const double * restrict voffluxx = interface->voffluxx;
  const double * restrict voffluxy = interface->voffluxy;
  double * restrict vofsrca = interface->vofsrca;
  double * restrict vofsrcb = interface->vofsrcb;
  memcpy(vofsrcb, vofsrca, VOFSRCA_SIZE_0 * VOFSRCA_SIZE_1 * sizeof(double));
  /* ! compute right-hand-side of advection equation ! 8 ! */
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      VOFSRCA(i, j) =
        -(-VOFFLUXX(i, j)+VOFFLUXX(i+1, j  ))/DXF(i)
        -(-VOFFLUXY(i, j)+VOFFLUXY(i  , j+1))/dy
      ;
    }
  }
  return 0;
}

static int interface_advect_vof(const domain_t *domain, const int rkstep, const double dt, interface_t *interface){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double alpha = RKCOEFS[rkstep].alpha;
  const double beta  = RKCOEFS[rkstep].beta;
  const double * restrict vofsrca = interface->vofsrca;
  const double * restrict vofsrcb = interface->vofsrcb;
  double * restrict vof = interface->vof;
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      VOF(i, j) += (
          + alpha * VOFSRCA(i, j)
          + beta  * VOFSRCB(i, j)
      )*dt;
    }
  }
  interface_update_boundaries_vof(domain, vof);
  return 0;
}


int interface_update_vof(const domain_t *domain, const int rkstep, const double dt, const fluid_t *fluid, interface_t *interface){
  interface_compute_flux_x(domain, fluid, interface);
  interface_compute_flux_y(domain, fluid, interface);
  interface_compute_rhs(domain, interface);
  interface_advect_vof(domain, rkstep, dt, interface);
  return 0;
}

