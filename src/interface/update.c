#include <string.h>
#include "common.h"
#include "domain.h"
#include "fluid.h"
#include "interface.h"
#include "arrays/fluid.h"
#include "arrays/interface.h"
#include "internal.h"


#if NDIMS == 2

static int interface_compute_flux_x(const domain_t *domain, const fluid_t *fluid, interface_t *interface){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict ux = fluid->ux;
  const double * restrict gps = interface->gps;
  const double * restrict gws = interface->gws;
  const double * restrict vof = interface->vof;
  const nrml_t * restrict normal = interface->normal;
  double * restrict voffluxx = interface->voffluxx;
  memset(voffluxx, 0, VOFFLUXX_SIZE_0 * VOFFLUXX_SIZE_1 * sizeof(double));
  for(int j = 1; j <= jsize; j++){
    for(int i = 2; i <= isize; i++){
      int ii;
      double x;
      if(UX(i, j) < 0.){
        ii = i;
        x = -0.5;
      }else{
        ii = i-1;
        x = +0.5;
      }
      if(VOF(ii, j) < VOFMIN || 1.-VOFMIN < VOF(ii, j)){
        VOFFLUXX(i, j) = VOF(ii, j);
      }else{
        double a10 = NORMAL(ii, j).a10;
        double a01 = NORMAL(ii, j).a01;
        double a00 = NORMAL(ii, j).a00;
        for(int ng = 0; ng < ORDER_GAUSS; ng++){
          double gw = gws[ng];
          double gy = gps[ng];
          VOFFLUXX(i, j) += gw*H(a00, a10, a01, x, gy);
        }
      }
      VOFFLUXX(i, j) *= UX(i, j);
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
  memset(voffluxy, 0, VOFFLUXY_SIZE_0 * VOFFLUXY_SIZE_1 * sizeof(double));
  for(int j = 1; j <= jsize+1; j++){
    for(int i = 1; i <= isize; i++){
      int jj;
      double y;
      if(UY(i, j) < 0.){
        jj = j;
        y = -0.5;
      }else{
        jj = j-1;
        y = +0.5;
      }
      if(VOF(i, jj) < VOFMIN || 1.-VOFMIN < VOF(i, jj)){
        VOFFLUXY(i, j) = VOF(i, jj);
      }else{
        double a10 = NORMAL(i, jj).a10;
        double a01 = NORMAL(i, jj).a01;
        double a00 = NORMAL(i, jj).a00;
        for(int ng = 0; ng < ORDER_GAUSS; ng++){
          double gw = gws[ng];
          double gx = gps[ng];
          VOFFLUXY(i, j) += gw*H(a00, a10, a01, gx, y);
        }
      }
      VOFFLUXY(i, j) *= UY(i, j);
    }
  }
  return 0;
}

static int interface_compute_rhs(const domain_t *domain, interface_t *interface){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double dx = domain->dx;
  const double dy = domain->dy;
  const double * restrict voffluxx = interface->voffluxx;
  const double * restrict voffluxy = interface->voffluxy;
  double * restrict vofsrca = interface->vofsrca;
  double * restrict vofsrcb = interface->vofsrcb;
  memcpy(vofsrcb, vofsrca, VOFSRCA_SIZE_0 * VOFSRCA_SIZE_1 * sizeof(double));
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      VOFSRCA(i, j) =
        -(-VOFFLUXX(i, j)+VOFFLUXX(i+1, j  ))/dx
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

#else // NDIMS == 3

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
        int ii;
        double x;
        if(UX(i, j, k) < 0.){
          ii = i;
          x = -0.5;
        }else{
          ii = i-1;
          x = +0.5;
        }
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
        int jj;
        double y;
        if(UY(i, j, k) < 0.){
          jj = j;
          y = -0.5;
        }else{
          jj = j-1;
          y = +0.5;
        }
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
        int kk;
        double z;
        if(UZ(i, j, k) < 0.){
          kk = k;
          z = -0.5;
        }else{
          kk = k-1;
          z = +0.5;
        }
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
  const double dx = domain->dx;
  const double dy = domain->dy;
  const double dz = domain->dz;
  const double * restrict voffluxx = interface->voffluxx;
  const double * restrict voffluxy = interface->voffluxy;
  const double * restrict voffluxz = interface->voffluxz;
  double * restrict vofsrca = interface->vofsrca;
  double * restrict vofsrcb = interface->vofsrcb;
  memcpy(vofsrcb, vofsrca, VOFSRCA_SIZE_0 * VOFSRCA_SIZE_1 * VOFSRCA_SIZE_2 * sizeof(double));
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        VOFSRCA(i, j, k) =
          -(-VOFFLUXX(i, j, k)+VOFFLUXX(i+1, j  , k  ))/dx
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

#endif // NDIMS

int interface_update_vof(const domain_t *domain, const int rkstep, const double dt, const fluid_t *fluid, interface_t *interface){
  interface_compute_flux_x(domain, fluid, interface);
  interface_compute_flux_y(domain, fluid, interface);
#if NDIMS == 3
  interface_compute_flux_z(domain, fluid, interface);
#endif
  interface_compute_rhs(domain, interface);
  interface_advect_vof(domain, rkstep, dt, interface);
  return 0;
}

