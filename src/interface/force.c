#include "config.h"
#include "domain.h"
#include "interface.h"
#include "arrays/interface.h"



int interface_compute_surface_force(const domain_t *domain, interface_t *interface){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict dxc = domain->dxc;
  const double            dy  = domain->dy;
  const double We = config.get_double("We");
  const double * restrict vof = interface->vof;
  const double * restrict curv = interface->curv;
  double * restrict voffrcx = interface->voffrcx;
  double * restrict voffrcy = interface->voffrcy;
  /* ! compute surface tension force in x direction ! 7 ! */
  for(int j = 1; j <= jsize; j++){
    for(int i = 2; i <= isize; i++){
      double grad = (-VOF(i-1, j)+VOF(i, j))/DXC(i);
      double kappa = 0.5*(CURV(i-1, j)+CURV(i, j));
      VOFFRCX(i, j) = 1./We*grad*kappa;
    }
  }
  /* ! compute surface tension force in y direction ! 7 ! */
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      double grad = (-VOF(i, j-1)+VOF(i, j))/dy;
      double kappa = 0.5*(CURV(i, j-1)+CURV(i, j));
      VOFFRCY(i, j) = 1./We*grad*kappa;
    }
  }
  return 0;
}

