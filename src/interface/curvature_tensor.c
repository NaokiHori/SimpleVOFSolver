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
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict dxc = domain->dxc;
  const double            dy  = domain->dy;
#if NDIMS == 3
  const double            dz  = domain->dz;
#endif
  const double * restrict vof = interface->vof.data;
  vector_t * restrict dvof = interface->dvof.data;
#if NDIMS == 2
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
#else
  for(int k = 0; k <= ksize + 2; k++){
    for(int j = 0; j <= jsize + 2; j++){
      for(int i = 1; i <= isize + 1; i++){
        const double dx = DXC(i  );
        // x gradient | 6
        const double dvofdx = 1. / dx * (
            - VOF(i-1, j-1, k-1) + VOF(i  , j-1, k-1)
            - VOF(i-1, j  , k-1) + VOF(i  , j  , k-1)
            - VOF(i-1, j-1, k  ) + VOF(i  , j-1, k  )
            - VOF(i-1, j  , k  ) + VOF(i  , j  , k  )
        );
        // y gradient | 6
        const double dvofdy = 1. / dy * (
            - VOF(i-1, j-1, k-1) - VOF(i  , j-1, k-1)
            + VOF(i-1, j  , k-1) + VOF(i  , j  , k-1)
            - VOF(i-1, j-1, k  ) - VOF(i  , j-1, k  )
            + VOF(i-1, j  , k  ) + VOF(i  , j  , k  )
        );
        // z gradient | 6
        const double dvofdz = 1. / dz * (
            - VOF(i-1, j-1, k-1) - VOF(i  , j-1, k-1)
            - VOF(i-1, j  , k-1) - VOF(i  , j  , k-1)
            + VOF(i-1, j-1, k  ) + VOF(i  , j-1, k  )
            + VOF(i-1, j  , k  ) + VOF(i  , j  , k  )
        );
        // normalise and obtain corner normals | 9
        const double norm = sqrt(
            + pow(dvofdx, 2.)
            + pow(dvofdy, 2.)
            + pow(dvofdz, 2.)
        );
        const double norminv = 1. / fmax(norm, DBL_EPSILON);
        DVOF(i, j, k)[0] = dvofdx * norminv;
        DVOF(i, j, k)[1] = dvofdy * norminv;
        DVOF(i, j, k)[2] = dvofdz * norminv;
      }
    }
  }
#endif
  return 0;
}

static double compute_intercept(
    const double vof,
    const double normal[NDIMS]
){
  // Newton-Raphson method, loop terminating conditions | 2
  const int cntmax = 8;
  const double resmax = 1.e-12;
#if NDIMS == 2
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
#else
  // compute constants a priori | 14
  double exps[NGAUSS * NGAUSS * NGAUSS] = {0.};
  for(int kk = 0; kk < NGAUSS; kk++){
    for(int jj = 0; jj < NGAUSS; jj++){
      for(int ii = 0; ii < NGAUSS; ii++){
        exps[(kk * NGAUSS + jj) * NGAUSS + ii] = exp(
            -2. * vofbeta * (
              + normal[0] * gauss_ps[ii]
              + normal[1] * gauss_ps[jj]
              + normal[2] * gauss_ps[kk]
            )
        );
      }
    }
  }
#endif
  // initial guess | 1
  double val = 1. / vof - 1.;
  for(int cnt = 0; cnt < cntmax; cnt++){
    // sum up | 25
    double f0 = -1. * vof;
    double f1 = 0.;
#if NDIMS == 2
    for(int jj = 0; jj < NGAUSS; jj++){
      for(int ii = 0; ii < NGAUSS; ii++){
        const double weight = gauss_ws[ii] * gauss_ws[jj];
        const double p = exps[jj * NGAUSS + ii];
        const double denom = 1. / (1. + p * val);
        f0 += weight     * denom;
        f1 -= weight * p * denom * denom;
      }
    }
#else
    for(int kk = 0; kk < NGAUSS; kk++){
      for(int jj = 0; jj < NGAUSS; jj++){
        for(int ii = 0; ii < NGAUSS; ii++){
          const double weight = gauss_ws[ii] * gauss_ws[jj] * gauss_ws[kk];
          const double p = exps[(kk * NGAUSS + jj) * NGAUSS + ii];
          const double denom = 1. / (1. + p * val);
          f0 += weight     * denom;
          f1 -= weight * p * denom * denom;
        }
      }
    }
#endif
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
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict dxf = domain->dxf;
  const double            dy  = domain->dy;
#if NDIMS == 3
  const double            dz  = domain->dz;
#endif
  const double * restrict vof = interface->vof.data;
  const vector_t * restrict dvof = interface->dvof.data;
  normal_t * restrict normal = interface->normal.data;
#if NDIMS == 2
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
#else
  for(int k = 0; k <= ksize + 1; k++){
    for(int j = 0; j <= jsize + 1; j++){
      for(int i = 1; i <= isize; i++){
        const double dx = DXF(i  );
        const double lvof = VOF(i, j, k);
        // for (almost) single-phase region,
        //   surface reconstruction is not needed
        if(lvof < vofmin || 1. - vofmin < lvof){
          continue;
        }
        // average nx | 6
        double nx = (
            + DVOF(i  , j  , k  )[0] + DVOF(i+1, j  , k  )[0]
            + DVOF(i  , j+1, k  )[0] + DVOF(i+1, j+1, k  )[0]
            + DVOF(i  , j  , k+1)[0] + DVOF(i+1, j  , k+1)[0]
            + DVOF(i  , j+1, k+1)[0] + DVOF(i+1, j+1, k+1)[0]
        );
        // average ny | 6
        double ny = (
            + DVOF(i  , j  , k  )[1] + DVOF(i+1, j  , k  )[1]
            + DVOF(i  , j+1, k  )[1] + DVOF(i+1, j+1, k  )[1]
            + DVOF(i  , j  , k+1)[1] + DVOF(i+1, j  , k+1)[1]
            + DVOF(i  , j+1, k+1)[1] + DVOF(i+1, j+1, k+1)[1]
        );
        // average nz | 6
        double nz = (
            + DVOF(i  , j  , k  )[2] + DVOF(i+1, j  , k  )[2]
            + DVOF(i  , j+1, k  )[2] + DVOF(i+1, j+1, k  )[2]
            + DVOF(i  , j  , k+1)[2] + DVOF(i+1, j  , k+1)[2]
            + DVOF(i  , j+1, k+1)[2] + DVOF(i+1, j+1, k+1)[2]
        );
        // normalise and obtain center normals | 12
        nx = nx / dx;
        ny = ny / dy;
        nz = nz / dz;
        const double norm = sqrt(
            + pow(nx, 2.)
            + pow(ny, 2.)
            + pow(nz, 2.)
        );
        const double norminv = 1. / fmax(norm, DBL_EPSILON);
        nx *= norminv;
        ny *= norminv;
        nz *= norminv;
        const double seg = compute_intercept(
            lvof,
            (const double [NDIMS]){nx, ny, nz}
        );
        // store normal and intercept | 4
        NORMAL(i, j, k)[0] = nx;
        NORMAL(i, j, k)[1] = ny;
        NORMAL(i, j, k)[2] = nz;
        NORMAL(i, j, k)[3] = seg;
      }
    }
  }
#endif
  return 0;
}

static int compute_curvature(
    const domain_t * domain,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict dxf = domain->dxf;
  const double            dy  = domain->dy;
#if NDIMS == 3
  const double            dz  = domain->dz;
#endif
  const vector_t * restrict dvof = interface->dvof.data;
  double * restrict curv = interface->curv.data;
#if NDIMS == 2
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
#else
  for(int k = 0; k <= ksize + 1; k++){
    for(int j = 0; j <= jsize + 1; j++){
      for(int i = 1; i <= isize; i++){
        const double dx = DXF(i  );
        // compute mean curvature from corner normals | 23
        const double dnxdx = 1. / dx * (
            - DVOF(i  , j  , k  )[0] + DVOF(i+1, j  , k  )[0]
            - DVOF(i  , j+1, k  )[0] + DVOF(i+1, j+1, k  )[0]
            - DVOF(i  , j  , k+1)[0] + DVOF(i+1, j  , k+1)[0]
            - DVOF(i  , j+1, k+1)[0] + DVOF(i+1, j+1, k+1)[0]
        );
        const double dnydy = 1. / dy * (
            - DVOF(i  , j  , k  )[1] - DVOF(i+1, j  , k  )[1]
            + DVOF(i  , j+1, k  )[1] + DVOF(i+1, j+1, k  )[1]
            - DVOF(i  , j  , k+1)[1] - DVOF(i+1, j  , k+1)[1]
            + DVOF(i  , j+1, k+1)[1] + DVOF(i+1, j+1, k+1)[1]
        );
        const double dnzdz = 1. / dz * (
            - DVOF(i  , j  , k  )[2] - DVOF(i+1, j  , k  )[2]
            - DVOF(i  , j+1, k  )[2] - DVOF(i+1, j+1, k  )[2]
            + DVOF(i  , j  , k+1)[2] + DVOF(i+1, j  , k+1)[2]
            + DVOF(i  , j+1, k+1)[2] + DVOF(i+1, j+1, k+1)[2]
        );
        CURV(i, j, k) = 0.25 * (
          - dnxdx
          - dnydy
          - dnzdz
        );
      }
    }
  }
#endif
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

