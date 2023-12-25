#include <math.h>
#include <float.h>
#include "param.h"
#include "domain.h"
#include "interface.h"
#include "internal.h"
#include "array_macros/domain/dxf.h"
#include "array_macros/domain/dxc.h"
#include "array_macros/interface/vof.h"
#include "array_macros/interface/dvof.h"
#include "array_macros/interface/normal.h"
#include "array_macros/interface/curv.h"

static int compute_gradient(
    const domain_t * domain,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict dxc = domain->dxc;
  const double            dy  = domain->dy;
  const double * restrict vof = interface->vof.data;
  vector_t * restrict dvof = interface->dvof.data;
  for(int j = 0; j <= jsize + 2; j++){
    for(int i = 1; i <= isize + 1; i++){
      // x gradient | 5
      const double dx = DXC(i  );
      const double dvofdx = 1. / dx * (
          - VOF(i-1, j-1) + VOF(i  , j-1)
          - VOF(i-1, j  ) + VOF(i  , j  )
      );
      // y gradient | 4
      const double dvofdy = 1. / dy * (
          - VOF(i-1, j-1) - VOF(i  , j-1)
          + VOF(i-1, j  ) + VOF(i  , j  )
      );
      // normalise and obtain corner normals | 7
      const double norm = sqrt(
          + pow(dvofdx, 2.)
          + pow(dvofdy, 2.)
      );
      const double norminv = 1. / fmax(norm, DBL_EPSILON);
      DVOF(i, j)[0] = dvofdx * norminv;
      DVOF(i, j)[1] = dvofdy * norminv;
    }
  }
  return 0;
}

static double compute_intercept(
    const double vof,
    const double normal[NDIMS]
){
  // Newton-Raphson method, loop terminating conditions | 2
  const int cntmax = 8;
  const double resmax = 1.e-12;
  // compute constants a priori | 11
  double exps[NGAUSS * NGAUSS] = {0.};
  for(int jj = 0; jj < NGAUSS; jj++){
    for(int ii = 0; ii < NGAUSS; ii++){
      exps[jj * NGAUSS + ii] = exp(
          -2. * vofbeta * (
            + normal[0] * gauss_ps[ii]
            + normal[1] * gauss_ps[jj]
          )
      );
    }
  }
  // initial guess | 1
  double val = 1. / vof - 1.;
  for(int cnt = 0; cnt < cntmax; cnt++){
    // sum up | 25
    double f0 = -1. * vof;
    double f1 = 0.;
    for(int jj = 0; jj < NGAUSS; jj++){
      for(int ii = 0; ii < NGAUSS; ii++){
        const double weight = gauss_ws[ii] * gauss_ws[jj];
        const double p = exps[jj * NGAUSS + ii];
        const double denom = 1. / (1. + p * val);
        f0 += weight     * denom;
        f1 -= weight * p * denom * denom;
      }
    }
    // update D | 1
    val -= f0 / f1;
    if(fabs(f0) < resmax){
      break;
    }
  }
  // convert D to d | 1
  return -0.5 / vofbeta * log(val);
}

static int compute_normal(
    const domain_t * domain,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict dxf = domain->dxf;
  const double            dy  = domain->dy;
  const double * restrict vof = interface->vof.data;
  const vector_t * restrict dvof = interface->dvof.data;
  normal_t * restrict normal = interface->normal.data;
  for(int j = 0; j <= jsize + 1; j++){
    for(int i = 1; i <= isize; i++){
      const double dx = DXF(i  );
      const double lvof = VOF(i, j);
      // for (almost) single-phase region,
      //   surface reconstruction is not needed
      if(lvof < vofmin || 1. - vofmin < lvof){
        continue;
      }
      // average nx | 4
      double nx = (
          + DVOF(i  , j  )[0] + DVOF(i+1, j  )[0]
          + DVOF(i  , j+1)[0] + DVOF(i+1, j+1)[0]
      );
      // average ny | 4
      double ny = (
          + DVOF(i  , j  )[1] + DVOF(i+1, j  )[1]
          + DVOF(i  , j+1)[1] + DVOF(i+1, j+1)[1]
      );
      // normalise and obtain center normals | 9
      nx /= dx;
      ny /= dy;
      const double norm = sqrt(
          + pow(nx, 2.)
          + pow(ny, 2.)
      );
      const double norminv = 1. / fmax(norm, DBL_EPSILON);
      nx *= norminv;
      ny *= norminv;
      const double seg = compute_intercept(
          lvof,
          (const double [NDIMS]){nx, ny}
      );
      // store normal and intercept | 3
      NORMAL(i, j)[0] = nx;
      NORMAL(i, j)[1] = ny;
      NORMAL(i, j)[2] = seg;
    }
  }
  return 0;
}

static int compute_curvature(
    const domain_t * domain,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict dxf = domain->dxf;
  const double            dy  = domain->dy;
  const vector_t * restrict dvof = interface->dvof.data;
  double * restrict curv = interface->curv.data;
  for(int j = 0; j <= jsize + 1; j++){
    for(int i = 1; i <= isize; i++){
      const double dx = DXF(i  );
      // compute mean curvature from corner normals | 12
      const double dnxdx = 1. / dx * (
          - DVOF(i  , j  )[0] + DVOF(i+1, j  )[0]
          - DVOF(i  , j+1)[0] + DVOF(i+1, j+1)[0]
      );
      const double dnydy = 1. / dy * (
          - DVOF(i  , j  )[1] - DVOF(i+1, j  )[1]
          + DVOF(i  , j+1)[1] + DVOF(i+1, j+1)[1]
      );
      CURV(i, j) = 0.5 * (
        - dnxdx
        - dnydy
      );
    }
  }
  return 0;
}

int interface_compute_curvature_tensor(
    const domain_t * domain,
    interface_t * interface
){
  compute_gradient(domain, interface);
  compute_normal(domain, interface);
  compute_curvature(domain, interface);
  return 0;
}

