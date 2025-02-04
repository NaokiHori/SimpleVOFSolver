#include "domain.h"
#include "interface.h"
#include "internal.h"
#include "array_macros/domain/dxc.h"
#include "array_macros/interface/vof.h"
#include "array_macros/interface/ifrcx.h"
#include "array_macros/interface/ifrcy.h"
#include "array_macros/interface/ifrcz.h"
#include "array_macros/interface/curv.h"

static int compute_force_x(
    const domain_t * domain,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict dxc = domain->dxc;
  const double tension = interface->tension;
  const double * restrict vof = interface->vof.data;
  const double * restrict curv = interface->curv.data;
  double * restrict ifrcx = interface->ifrcx.data;
  for(int j = 1; j <= jsize; j++){
    for(int i = 2; i <= isize; i++){
      // compute surface tension force in x direction | 10
      const double dx = DXC(i  );
      const double grad = 1. / dx * (
          - VOF(i-1, j  )
          + VOF(i  , j  )
      );
      const double kappa = 0.5 * (
          + CURV(i-1, j  )
          + CURV(i  , j  )
      );
      IFRCX(i, j) = tension * grad * kappa;
    }
  }
  return 0;
}

static int compute_force_y(
    const domain_t * domain,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double dy = domain->dy;
  const double tension = interface->tension;
  const double * restrict vof = interface->vof.data;
  const double * restrict curv = interface->curv.data;
  double * restrict ifrcy = interface->ifrcy.data;
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      // compute surface tension force in y direction | 9
      const double grad = 1. / dy * (
          - VOF(i  , j-1)
          + VOF(i  , j  )
      );
      const double kappa = 0.5 * (
          + CURV(i  , j-1)
          + CURV(i  , j  )
      );
      IFRCY(i, j) = tension * grad * kappa;
    }
  }
  return 0;
}


int interface_compute_force(
    const domain_t * domain,
    interface_t * interface
){
  compute_force_x(domain, interface);
  compute_force_y(domain, interface);
  return 0;
}

