#if NDIMS == 3
#include "domain.h"
#include "fluid.h"
#include "interface.h"
#include "../internal.h"
#include "internal.h"
#include "array_macros/fluid/uz.h"
#include "array_macros/interface/vof.h"
#include "array_macros/interface/normal.h"
#include "array_macros/interface/flxz.h"

int compute_flux_z(
    const domain_t * domain,
    const fluid_t * fluid,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict uz = fluid->uz.data;
  const double * restrict vof = interface->vof.data;
  const normal_t * restrict normal = interface->normal.data;
  double * restrict flxz = interface->flxz.data;
  for(int k = 1; k <= ksize + 1; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // use upwind information | 3
        const double vel = UZ(i, j, k);
        const int    kk = vel < 0. ?    k : k - 1;
        const double  z = vel < 0. ? -0.5 :  +0.5;
        // evaluate flux | 15
        const double lvof = VOF(i, j, kk);
        if(lvof < vofmin || 1. - vofmin < lvof){
          FLXZ(i, j, k) = vel * lvof;
          continue;
        }
        double flux = 0.;
        for(int jj = 0; jj < NGAUSS; jj++){
          for(int ii = 0; ii < NGAUSS; ii++){
            const double w = gauss_ws[ii] * gauss_ws[jj];
            const double x = gauss_ps[ii];
            const double y = gauss_ps[jj];
            flux += w * indicator(NORMAL(i, j, kk), (const double [NDIMS]){x, y, z});
          }
        }
        FLXZ(i, j, k) = vel * flux;
      }
    }
  }
  return 0;
}
#endif
