#include "domain.h"
#include "fluid.h"
#include "interface.h"
#include "../internal.h"
#include "internal.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/interface/vof.h"
#include "array_macros/interface/normal.h"
#include "array_macros/interface/flxy.h"

int compute_flux_y(
    const domain_t * domain,
    const fluid_t * fluid,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict uy = fluid->uy.data;
  const double * restrict vof = interface->vof.data;
  const normal_t * restrict normal = interface->normal.data;
  double * restrict flxy = interface->flxy.data;
  for(int j = 1; j <= jsize + 1; j++){
    for(int i = 1; i <= isize; i++){
      // use upwind information | 3
      const double vel = UY(i, j);
      const int    jj = vel < 0. ?    j : j - 1;
      const double  y = vel < 0. ? -0.5 :  +0.5;
      // evaluate flux | 12
      const double lvof = VOF(i, jj);
      if(lvof < vofmin || 1. - vofmin < lvof){
        FLXY(i, j) = vel * lvof;
        continue;
      }
      double flux = 0.;
      for(int ii = 0; ii < NGAUSS; ii++){
        const double w = gauss_ws[ii];
        const double x = gauss_ps[ii];
        flux += w * indicator(NORMAL(i, jj), (const double [NDIMS]){x, y});
      }
      FLXY(i, j) = vel * flux;
    }
  }
  return 0;
}

