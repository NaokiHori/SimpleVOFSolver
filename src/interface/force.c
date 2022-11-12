#include "config.h"
#include "domain.h"
#include "interface.h"
#include "arrays/interface.h"


#if NDIMS == 2

int interface_compute_surface_force(const domain_t *domain, interface_t *interface){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double dx = domain->dx;
  const double dy = domain->dy;
  const double We = config.get_double("We");
  const double * restrict vof = interface->vof;
  const double * restrict curv = interface->curv;
  double * restrict voffrcx = interface->voffrcx;
  double * restrict voffrcy = interface->voffrcy;
  for(int j = 1; j <= jsize; j++){
    for(int i = 2; i <= isize; i++){
      double grad = (-VOF(i-1, j)+VOF(i, j))/dx;
      double kappa = 0.5*(CURV(i-1, j)+CURV(i, j));
      VOFFRCX(i, j) = 1./We*grad*kappa;
    }
  }
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      double grad = (-VOF(i, j-1)+VOF(i, j))/dy;
      double kappa = 0.5*(CURV(i, j-1)+CURV(i, j));
      VOFFRCY(i, j) = 1./We*grad*kappa;
    }
  }
  return 0;
}

#else

int interface_compute_surface_force(const domain_t *domain, interface_t *interface){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double dx = domain->dx;
  const double dy = domain->dy;
  const double dz = domain->dz;
  const double We = config.get_double("We");
  const double * restrict vof = interface->vof;
  const double * restrict curv = interface->curv;
  double * restrict voffrcx = interface->voffrcx;
  double * restrict voffrcy = interface->voffrcy;
  double * restrict voffrcz = interface->voffrcz;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        double grad = (-VOF(i-1, j, k)+VOF(i, j, k))/dx;
        double kappa = 0.5*(CURV(i-1, j, k)+CURV(i, j, k));
        VOFFRCX(i, j, k) = 1./We*grad*kappa;
      }
    }
  }
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        double grad = (-VOF(i, j-1, k)+VOF(i, j, k))/dy;
        double kappa = 0.5*(CURV(i, j-1, k)+CURV(i, j, k));
        VOFFRCY(i, j, k) = 1./We*grad*kappa;
      }
    }
  }
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        double grad = (-VOF(i, j, k-1)+VOF(i, j, k))/dz;
        double kappa = 0.5*(CURV(i, j, k-1)+CURV(i, j, k));
        VOFFRCZ(i, j, k) = 1./We*grad*kappa;
      }
    }
  }
  return 0;
}

#endif
