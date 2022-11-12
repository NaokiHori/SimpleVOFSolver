#include <math.h>
#include <float.h>
#include "domain.h"
#include "interface.h"
#include "internal.h"
#include "arrays/interface.h"


#if NDIMS == 2

static int compute_gradient(const domain_t *domain, interface_t *interface){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict vof = interface->vof;
  dvof_t * restrict dvof = interface->dvof;
  for(int j = 0; j <= jsize+2; j++){
    for(int i = 1; i <= isize+1; i++){
      // x gradient
      double dxi = i == 1 || i == isize+1 ? 2. : 1.;
      double dvofdx = 0.5*dxi*(
          -VOF(i-1, j-1)+VOF(i  , j-1)
          -VOF(i-1, j  )+VOF(i  , j  )
      );
      // y gradient
      double dyi = 1.;
      double dvofdy = 0.5*dyi*(
          -VOF(i-1, j-1)-VOF(i  , j-1)
          +VOF(i-1, j  )+VOF(i  , j  )
      );
      // normalise
      double norm = sqrt(
          +pow(dvofdx, 2.)
          +pow(dvofdy, 2.)
      );
      norm = fmax(norm, DBL_EPSILON);
      DVOF(i, j).x = dvofdx/norm;
      DVOF(i, j).y = dvofdy/norm;
    }
  }
  return 0;
}

static int compute_normal(const domain_t *domain, interface_t *interface){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict vof = interface->vof;
  const dvof_t * restrict dvof = interface->dvof;
  nrml_t * restrict normal = interface->normal;
  for(int j = 0; j <= jsize+1; j++){
    for(int i = 1; i <= isize; i++){
      double nx = 0.25*(
          +DVOF(i  , j  ).x
          +DVOF(i+1, j  ).x
          +DVOF(i  , j+1).x
          +DVOF(i+1, j+1).x
      );
      double ny = 0.25*(
          +DVOF(i  , j  ).y
          +DVOF(i+1, j  ).y
          +DVOF(i  , j+1).y
          +DVOF(i+1, j+1).y
      );
      double norm = sqrt(
          +pow(nx, 2.)
          +pow(ny, 2.)
      );
      norm = fmax(norm, DBL_EPSILON);
      nx /= norm;
      ny /= norm;
      NORMAL(i, j).a10 = nx;
      NORMAL(i, j).a01 = ny;
      // Newton-Raphson method to obtain distance (nrans)
      {
        // max iteration
        const int nrcntmax = 10;
        // max residual
        const double nreps = 1.e-15;
        // used Gaussian quadratures
        const double * restrict gps = interface->gps;
        const double * restrict gws = interface->gws;
        double exp_psis[ORDER_GAUSS][ORDER_GAUSS];
        for(int jj = 0; jj < ORDER_GAUSS; jj++){
          for(int ii = 0; ii < ORDER_GAUSS; ii++){
            exp_psis[jj][ii] = exp(-2.*VOFBETA*psi(nx, ny, gps[ii], gps[jj]));
          }
        }
        // solution of iteration
        double nrans = 1. / VOF(i, j) - 1.;
        int cnt = 0;
        double f, fp;
        do {
          f = -VOF(i, j);
          fp = 0.;
          for(int jj = 0; jj < ORDER_GAUSS; jj++){
            for(int ii = 0; ii < ORDER_GAUSS; ii++){
              double tmp = gws[jj]*gws[ii];
              // function
              f  += tmp/(1.+exp_psis[jj][ii]*nrans);
              // its derivative (f^{prime})
              fp -= tmp*exp_psis[jj][ii]/pow(1.+exp_psis[jj][ii]*nrans, 2.);
            }
          }
          nrans -= f/fp;
          cnt++;
          if(cnt > nrcntmax){
            break;
          }
        } while(fabs(f) > nreps);
        NORMAL(i, j).a00 = -0.5/VOFBETA*log(nrans);
      }
    }
  }
  return 0;
}

static int compute_curvature(const domain_t *domain, interface_t *interface){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double dx = domain->dx;
  const double dy = domain->dy;
  const dvof_t * restrict dvof = interface->dvof;
  double * restrict curv = interface->curv;
  for(int j = 0; j <= jsize+1; j++){
    for(int i = 1; i <= isize; i++){
      double dnxdx = 0.5/dx*(
          - DVOF(i  , j  ).x + DVOF(i+1, j  ).x
          - DVOF(i  , j+1).x + DVOF(i+1, j+1).x
      );
      double dnydy = 0.5/dy*(
          - DVOF(i  , j  ).y - DVOF(i+1, j  ).y
          + DVOF(i  , j+1).y + DVOF(i+1, j+1).y
      );
      CURV(i, j) = -dnxdx-dnydy;
    }
  }
  return 0;
}

#else // NDIMS == 3

static int compute_gradient(const domain_t *domain, interface_t *interface){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict vof = interface->vof;
  dvof_t * restrict dvof = interface->dvof;
  for(int k = 0; k <= ksize+2; k++){
    for(int j = 0; j <= jsize+2; j++){
      for(int i = 1; i <= isize+1; i++){
        // x gradient
        double dxi = i == 1 || i == isize+1 ? 2. : 1.;
        double dvofdx = 0.25*dxi*(
            -VOF(i-1, j-1, k-1)+VOF(i  , j-1, k-1)
            -VOF(i-1, j  , k-1)+VOF(i  , j  , k-1)
            -VOF(i-1, j-1, k  )+VOF(i  , j-1, k  )
            -VOF(i-1, j  , k  )+VOF(i  , j  , k  )
        );
        // y gradient
        double dyi = 1.;
        double dvofdy = 0.25*dyi*(
            -VOF(i-1, j-1, k-1)-VOF(i  , j-1, k-1)
            +VOF(i-1, j  , k-1)+VOF(i  , j  , k-1)
            -VOF(i-1, j-1, k  )-VOF(i  , j-1, k  )
            +VOF(i-1, j  , k  )+VOF(i  , j  , k  )
        );
        // z gradient
        double dzi = 1.;
        double dvofdz = 0.25*dzi*(
            -VOF(i-1, j-1, k-1)-VOF(i  , j-1, k-1)
            -VOF(i-1, j  , k-1)-VOF(i  , j  , k-1)
            +VOF(i-1, j-1, k  )+VOF(i  , j-1, k  )
            +VOF(i-1, j  , k  )+VOF(i  , j  , k  )
        );
        // normalise
        double norm = sqrt(
            +pow(dvofdx, 2.)
            +pow(dvofdy, 2.)
            +pow(dvofdz, 2.)
        );
        norm = fmax(norm, DBL_EPSILON);
        DVOF(i, j, k).x = dvofdx/norm;
        DVOF(i, j, k).y = dvofdy/norm;
        DVOF(i, j, k).z = dvofdz/norm;
      }
    }
  }
  return 0;
}

static int compute_normal(const domain_t *domain, interface_t *interface){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict vof = interface->vof;
  const dvof_t * restrict dvof = interface->dvof;
  nrml_t * restrict normal = interface->normal;
  for(int k = 0; k <= ksize+1; k++){
    for(int j = 0; j <= jsize+1; j++){
      for(int i = 1; i <= isize; i++){
        double nx = 0.125*(
            +DVOF(i  , j  , k  ).x
            +DVOF(i+1, j  , k  ).x
            +DVOF(i  , j+1, k  ).x
            +DVOF(i+1, j+1, k  ).x
            +DVOF(i  , j  , k+1).x
            +DVOF(i+1, j  , k+1).x
            +DVOF(i  , j+1, k+1).x
            +DVOF(i+1, j+1, k+1).x
        );
        double ny = 0.125*(
            +DVOF(i  , j  , k  ).y
            +DVOF(i+1, j  , k  ).y
            +DVOF(i  , j+1, k  ).y
            +DVOF(i+1, j+1, k  ).y
            +DVOF(i  , j  , k+1).y
            +DVOF(i+1, j  , k+1).y
            +DVOF(i  , j+1, k+1).y
            +DVOF(i+1, j+1, k+1).y
        );
        double nz = 0.125*(
            +DVOF(i  , j  , k  ).z
            +DVOF(i+1, j  , k  ).z
            +DVOF(i  , j+1, k  ).z
            +DVOF(i+1, j+1, k  ).z
            +DVOF(i  , j  , k+1).z
            +DVOF(i+1, j  , k+1).z
            +DVOF(i  , j+1, k+1).z
            +DVOF(i+1, j+1, k+1).z
        );
        double norm = sqrt(
            +pow(nx, 2.)
            +pow(ny, 2.)
            +pow(nz, 2.)
        );
        norm = fmax(norm, DBL_EPSILON);
        nx /= norm;
        ny /= norm;
        nz /= norm;
        NORMAL(i, j, k).a100 = nx;
        NORMAL(i, j, k).a010 = ny;
        NORMAL(i, j, k).a001 = nz;
        // Newton-Raphson method to obtain distance (nrans)
        {
          // max iteration
          const int nrcntmax = 10;
          // max residual
          const double nreps = 1.e-15;
          // used Gaussian quadratures
          const double * restrict gps = interface->gps;
          const double * restrict gws = interface->gws;
          double exp_psis[ORDER_GAUSS][ORDER_GAUSS][ORDER_GAUSS];
          for(int kk = 0; kk < ORDER_GAUSS; kk++){
            for(int jj = 0; jj < ORDER_GAUSS; jj++){
              for(int ii = 0; ii < ORDER_GAUSS; ii++){
                exp_psis[kk][jj][ii] = exp(-2.*VOFBETA*psi(nx, ny, nz, gps[ii], gps[jj], gps[kk]));
              }
            }
          }
          // solution of iteration
          double nrans = 1. / VOF(i, j, k) - 1.;
          int cnt = 0;
          double f, fp;
          do {
            f = -VOF(i, j, k);
            fp = 0.;
            for(int kk = 0; kk < ORDER_GAUSS; kk++){
              for(int jj = 0; jj < ORDER_GAUSS; jj++){
                for(int ii = 0; ii < ORDER_GAUSS; ii++){
                  double tmp = gws[ii] * gws[jj] * gws[kk];
                  // function
                  f  += tmp/(1.+exp_psis[kk][jj][ii]*nrans);
                  // its derivative (f^{prime})
                  fp -= tmp*exp_psis[kk][jj][ii]/pow(1.+exp_psis[kk][jj][ii]*nrans, 2.);
                }
              }
            }
            nrans -= f/fp;
            cnt++;
            if(cnt > nrcntmax){
              break;
            }
          } while(fabs(f) > nreps);
          NORMAL(i, j, k).a000 = -0.5/VOFBETA*log(nrans);
        }
      }
    }
  }
  return 0;
}

static int compute_curvature(const domain_t *domain, interface_t *interface){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double dx = domain->dx;
  const double dy = domain->dy;
  const double dz = domain->dz;
  const dvof_t * restrict dvof = interface->dvof;
  double * restrict curv = interface->curv;
  for(int k = 0; k <= ksize+1; k++){
    for(int j = 0; j <= jsize+1; j++){
      for(int i = 1; i <= isize; i++){
        double dnxdx = 0.25/dx*(
            - DVOF(i  , j  , k  ).x + DVOF(i+1, j  , k  ).x
            - DVOF(i  , j+1, k  ).x + DVOF(i+1, j+1, k  ).x
            - DVOF(i  , j  , k+1).x + DVOF(i+1, j  , k+1).x
            - DVOF(i  , j+1, k+1).x + DVOF(i+1, j+1, k+1).x
        );
        double dnydy = 0.25/dy*(
            - DVOF(i  , j  , k  ).y - DVOF(i+1, j  , k  ).y
            + DVOF(i  , j+1, k  ).y + DVOF(i+1, j+1, k  ).y
            - DVOF(i  , j  , k+1).y - DVOF(i+1, j  , k+1).y
            + DVOF(i  , j+1, k+1).y + DVOF(i+1, j+1, k+1).y
        );
        double dnzdz = 0.25/dz*(
            - DVOF(i  , j  , k  ).z - DVOF(i+1, j  , k  ).z
            - DVOF(i  , j+1, k  ).z - DVOF(i+1, j+1, k  ).z
            + DVOF(i  , j  , k+1).z + DVOF(i+1, j  , k+1).z
            + DVOF(i  , j+1, k+1).z + DVOF(i+1, j+1, k+1).z
        );
        CURV(i, j, k) = -dnxdx-dnydy-dnzdz;
      }
    }
  }
  return 0;
}

#endif // NDIMS

int interface_compute_curvature_tensor(const domain_t *domain, interface_t *interface){
  compute_gradient(domain, interface);
  compute_normal(domain, interface);
  compute_curvature(domain, interface);
  return 0;
}

