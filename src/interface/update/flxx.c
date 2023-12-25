#include "domain.h"
#include "fluid.h"
#include "interface.h"
#include "../internal.h"
#include "internal.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/interface/vof.h"
#include "array_macros/interface/normal.h"
#include "array_macros/interface/flxx.h"

int compute_flux_x(
    const domain_t * domain,
    const fluid_t * fluid,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict ux = fluid->ux.data;
  const double * restrict vof = interface->vof.data;
  const normal_t * restrict normal = interface->normal.data;
  double * restrict flxx = interface->flxx.data;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        // use upwind information | 3
        const double vel = UX(i, j, k);
        const int    ii = vel < 0. ?    i : i - 1;
        const double  x = vel < 0. ? -0.5 :  +0.5;
        // evaluate flux | 15
        const double lvof = VOF(ii, j, k);
        if(lvof < vofmin || 1. - vofmin < lvof){
          FLXX(i, j, k) = vel * lvof;
          continue;
        }
        double flux = 0.;
        for(int kk = 0; kk < NGAUSS; kk++){
          for(int jj = 0; jj < NGAUSS; jj++){
            const double w = gauss_ws[jj] * gauss_ws[kk];
            const double y = gauss_ps[jj];
            const double z = gauss_ps[kk];
            flux += w * indicator(NORMAL(ii, j, k), (const double [NDIMS]){x, y, z});
          }
        }
        FLXX(i, j, k) = vel * flux;
      }
    }
  }
  return 0;
}

