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
  const int ksize = domain->mysizes[2];
  const double * restrict ux = fluid->ux;
  const double * restrict gps = interface->gps;
  const double * restrict gws = interface->gws;
  const double * restrict vof = interface->vof;
  const nrml_t * restrict normal = interface->normal;
  double * restrict voffluxx = interface->voffluxx;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        /* ! use upwind information ! 9 ! */
        int ii;
        double x;
        if(UX(i, j, k) < 0.){
          ii = i;
          x = -0.5;
        }else{
          ii = i-1;
          x = +0.5;
        }
        /* ! evaluate flux ! 15 ! */
        double flux = 0.;
        if(VOF(ii, j, k) < VOFMIN || 1.-VOFMIN < VOF(ii, j, k)){
          flux = VOF(ii, j, k);
        }else{
          double a100 = NORMAL(ii, j, k).a100;
          double a010 = NORMAL(ii, j, k).a010;
          double a001 = NORMAL(ii, j, k).a001;
          double a000 = NORMAL(ii, j, k).a000;
          for(int kk = 0; kk < ORDER_GAUSS; kk++){
            for(int jj = 0; jj < ORDER_GAUSS; jj++){
              flux += gws[jj] * gws[kk] * H(a000, a100, a010, a001, x, gps[jj], gps[kk]);
            }
          }
        }
        VOFFLUXX(i, j, k) = flux * UX(i, j, k);
      }
    }
  }
  return 0;
}

static int interface_compute_flux_y(const domain_t *domain, const fluid_t *fluid, interface_t *interface){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict uy = fluid->uy;
  const double * restrict gps = interface->gps;
  const double * restrict gws = interface->gws;
  const double * restrict vof = interface->vof;
  const nrml_t * restrict normal = interface->normal;
  double * restrict voffluxy = interface->voffluxy;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize+1; j++){
      for(int i = 1; i <= isize; i++){
        /* ! use upwind information ! 9 ! */
        int jj;
        double y;
        if(UY(i, j, k) < 0.){
          jj = j;
          y = -0.5;
        }else{
          jj = j-1;
          y = +0.5;
        }
        /* ! evaluate flux ! 15 ! */
        double flux = 0.;
        if(VOF(i, jj, k) < VOFMIN || 1.-VOFMIN < VOF(i, jj, k)){
          flux = VOF(i, jj, k);
        }else{
          double a100 = NORMAL(i, jj, k).a100;
          double a010 = NORMAL(i, jj, k).a010;
          double a001 = NORMAL(i, jj, k).a001;
          double a000 = NORMAL(i, jj, k).a000;
          for(int kk = 0; kk < ORDER_GAUSS; kk++){
            for(int ii = 0; ii < ORDER_GAUSS; ii++){
              flux += gws[ii] * gws[kk] * H(a000, a100, a010, a001, gps[ii], y, gps[kk]);
            }
          }
        }
        VOFFLUXY(i, j, k) = flux * UY(i, j, k);
      }
    }
  }
  return 0;
}

static int interface_compute_flux_z(const domain_t *domain, const fluid_t *fluid, interface_t *interface){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict uz = fluid->uz;
  const double * restrict gps = interface->gps;
  const double * restrict gws = interface->gws;
  const double * restrict vof = interface->vof;
  const nrml_t * restrict normal = interface->normal;
  double * restrict voffluxz = interface->voffluxz;
  for(int k = 1; k <= ksize+1; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        /* ! use upwind information ! 9 ! */
        int kk;
        double z;
        if(UZ(i, j, k) < 0.){
          kk = k;
          z = -0.5;
        }else{
          kk = k-1;
          z = +0.5;
        }
        /* ! evaluate flux ! 15 ! */
        double flux = 0.;
        if(VOF(i, j, kk) < VOFMIN || 1.-VOFMIN < VOF(i, j, kk)){
          flux = VOF(i, j, kk);
        }else{
          double a100 = NORMAL(i, j, kk).a100;
          double a010 = NORMAL(i, j, kk).a010;
          double a001 = NORMAL(i, j, kk).a001;
          double a000 = NORMAL(i, j, kk).a000;
          for(int jj = 0; jj < ORDER_GAUSS; jj++){
            for(int ii = 0; ii < ORDER_GAUSS; ii++){
              flux += gws[ii] * gws[jj] * H(a000, a100, a010, a001, gps[ii], gps[jj], z);
            }
          }
        }
        VOFFLUXZ(i, j, k) = flux * UZ(i, j, k);
      }
    }
  }
  return 0;
}

static int interface_compute_rhs(const domain_t *domain, interface_t *interface){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict dxf = domain->dxf;
  const double            dy  = domain->dy;
  const double            dz  = domain->dz;
  const double * restrict voffluxx = interface->voffluxx;
  const double * restrict voffluxy = interface->voffluxy;
  const double * restrict voffluxz = interface->voffluxz;
  double * restrict vofsrca = interface->vofsrca;
  double * restrict vofsrcb = interface->vofsrcb;
  memcpy(vofsrcb, vofsrca, VOFSRCA_SIZE_0 * VOFSRCA_SIZE_1 * VOFSRCA_SIZE_2 * sizeof(double));
  /* ! compute right-hand-side of advection equation ! 11 ! */
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        VOFSRCA(i, j, k) =
          -(-VOFFLUXX(i, j, k)+VOFFLUXX(i+1, j  , k  ))/DXF(i)
          -(-VOFFLUXY(i, j, k)+VOFFLUXY(i  , j+1, k  ))/dy
          -(-VOFFLUXZ(i, j, k)+VOFFLUXZ(i  , j  , k+1))/dz
        ;
      }
    }
  }
  return 0;
}

static int interface_advect_vof(const domain_t *domain, const int rkstep, const double dt, interface_t *interface){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double alpha = RKCOEFS[rkstep].alpha;
  const double beta  = RKCOEFS[rkstep].beta;
  const double * restrict vofsrca = interface->vofsrca;
  const double * restrict vofsrcb = interface->vofsrcb;
  double * restrict vof = interface->vof;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        VOF(i, j, k) += (
            + alpha * VOFSRCA(i, j, k)
            + beta  * VOFSRCB(i, j, k)
        )*dt;
      }
    }
  }
  interface_update_boundaries_vof(domain, vof);
  return 0;
}


int interface_update_vof(const domain_t *domain, const int rkstep, const double dt, const fluid_t *fluid, interface_t *interface){
  interface_compute_flux_x(domain, fluid, interface);
  interface_compute_flux_y(domain, fluid, interface);
  interface_compute_flux_z(domain, fluid, interface);
  interface_compute_rhs(domain, interface);
  interface_advect_vof(domain, rkstep, dt, interface);
  return 0;
}

