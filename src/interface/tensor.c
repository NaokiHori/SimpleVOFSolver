#include <math.h>
#include <float.h>
#include "domain.h"
#include "interface.h"
#include "internal.h"
#include "arrays/interface.h"


/* ! Newton-Raphson method, loop terminating conditions ! 4 ! */
// max iteration
static const int nrcntmax = 10;
// max residual
static const double nreps = 1.e-15;

#if NDIMS == 2

static int compute_gradient(const domain_t *domain, interface_t *interface){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict dxc = domain->dxc;
  const double            dy  = domain->dy;
  const double * restrict vof = interface->vof;
  dvof_t * restrict dvof = interface->dvof;
  for(int j = 0; j <= jsize+2; j++){
    for(int i = 1; i <= isize+1; i++){
      /* ! x gradient ! 4 ! */
      double dvofdx = 0.5/DXC(i)*(
          -VOF(i-1, j-1)+VOF(i  , j-1)
          -VOF(i-1, j  )+VOF(i  , j  )
      );
      /* ! y gradient ! 4 ! */
      double dvofdy = 0.5/dy*(
          -VOF(i-1, j-1)-VOF(i  , j-1)
          +VOF(i-1, j  )+VOF(i  , j  )
      );
      /* ! normalise and obtain corner normals ! 7 ! */
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
  const double * restrict dxf = domain->dxf;
  const double            dy  = domain->dy;
  const double * restrict vof = interface->vof;
  const dvof_t * restrict dvof = interface->dvof;
  nrml_t * restrict normal = interface->normal;
  // Gauss points
  const double * restrict gps = interface->gps;
  // Gauss weights
  const double * restrict gws = interface->gws;
  for(int j = 0; j <= jsize+1; j++){
    for(int i = 1; i <= isize; i++){
      // for (almost) single-phase region,
      //   surface reconstruction is not needed
      if(VOF(i, j) < VOFMIN || 1.-VOFMIN < VOF(i, j)){
        continue;
      }
      /* ! average nx ! 6 ! */
      double nx = 0.25*(
          +DVOF(i  , j  ).x
          +DVOF(i+1, j  ).x
          +DVOF(i  , j+1).x
          +DVOF(i+1, j+1).x
      );
      /* ! average ny ! 6 ! */
      double ny = 0.25*(
          +DVOF(i  , j  ).y
          +DVOF(i+1, j  ).y
          +DVOF(i  , j+1).y
          +DVOF(i+1, j+1).y
      );
      /* ! normalise and obtain center normals ! 9 ! */
      nx = nx / DXF(i);
      ny = ny / dy;
      double norm = sqrt(
          +pow(nx, 2.)
          +pow(ny, 2.)
      );
      norm = fmax(norm, DBL_EPSILON);
      nx /= norm;
      ny /= norm;
      /* Newton-Raphson iteration */
      double a00;
      {
        /* ! compute constants a priori ! 6 ! */
        double exp_psis[ORDER_GAUSS][ORDER_GAUSS];
        for(int jj = 0; jj < ORDER_GAUSS; jj++){
          for(int ii = 0; ii < ORDER_GAUSS; ii++){
            exp_psis[jj][ii] = exp(-2.*VOFBETA*psi(nx, ny, gps[ii], gps[jj]));
          }
        }
        /* ! initial guess ! 1 ! */
        a00 = 1. / VOF(i, j) - 1.;
        int cnt = 0;
        double f, fp;
        do {
          /* ! sum up ! 11 ! */
          f = -VOF(i, j);
          fp = 0.;
          for(int jj = 0; jj < ORDER_GAUSS; jj++){
            for(int ii = 0; ii < ORDER_GAUSS; ii++){
              double tmp = gws[jj]*gws[ii];
              // function
              f  += tmp/(1.+exp_psis[jj][ii]*a00);
              // its derivative (f^{prime})
              fp -= tmp*exp_psis[jj][ii]/pow(1.+exp_psis[jj][ii]*a00, 2.);
            }
          }
          /* ! update D ! 1 ! */
          a00 -= f/fp;
          cnt++;
          if(cnt > nrcntmax){
            break;
          }
        } while(fabs(f) > nreps);
        /* ! convert D to d ! 1 ! */
        a00 = -1./(2.*VOFBETA)*log(a00);
      }
      /* ! store normal and intercept ! 3 ! */
      NORMAL(i, j).a10 = nx;
      NORMAL(i, j).a01 = ny;
      NORMAL(i, j).a00 = a00;
    }
  }
  return 0;
}

static int compute_curvature(const domain_t *domain, interface_t *interface){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict dxf = domain->dxf;
  const double            dy  = domain->dy;
  const dvof_t * restrict dvof = interface->dvof;
  double * restrict curv = interface->curv;
  for(int j = 0; j <= jsize+1; j++){
    for(int i = 1; i <= isize; i++){
      /* ! compute mean curvature from corner normals ! 9 ! */
      double dnxdx = 0.5/DXF(i)*(
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
  const double * restrict dxc = domain->dxc;
  const double            dy  = domain->dy;
  const double            dz  = domain->dz;
  const double * restrict vof = interface->vof;
  dvof_t * restrict dvof = interface->dvof;
  for(int k = 0; k <= ksize+2; k++){
    for(int j = 0; j <= jsize+2; j++){
      for(int i = 1; i <= isize+1; i++){
        /* ! x gradient ! 6 ! */
        double dvofdx = 0.25/DXC(i)*(
            -VOF(i-1, j-1, k-1)+VOF(i  , j-1, k-1)
            -VOF(i-1, j  , k-1)+VOF(i  , j  , k-1)
            -VOF(i-1, j-1, k  )+VOF(i  , j-1, k  )
            -VOF(i-1, j  , k  )+VOF(i  , j  , k  )
        );
        /* ! y gradient ! 6 ! */
        double dvofdy = 0.25/dy*(
            -VOF(i-1, j-1, k-1)-VOF(i  , j-1, k-1)
            +VOF(i-1, j  , k-1)+VOF(i  , j  , k-1)
            -VOF(i-1, j-1, k  )-VOF(i  , j-1, k  )
            +VOF(i-1, j  , k  )+VOF(i  , j  , k  )
        );
        /* ! z gradient ! 6 ! */
        double dvofdz = 0.25/dz*(
            -VOF(i-1, j-1, k-1)-VOF(i  , j-1, k-1)
            -VOF(i-1, j  , k-1)-VOF(i  , j  , k-1)
            +VOF(i-1, j-1, k  )+VOF(i  , j-1, k  )
            +VOF(i-1, j  , k  )+VOF(i  , j  , k  )
        );
        /* ! normalise and obtain corner normals ! 9 ! */
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
  const double * restrict dxf = domain->dxf;
  const double            dy  = domain->dy;
  const double            dz  = domain->dz;
  const double * restrict vof = interface->vof;
  const dvof_t * restrict dvof = interface->dvof;
  nrml_t * restrict normal = interface->normal;
  // Gauss points
  const double * restrict gps = interface->gps;
  // Gauss weights
  const double * restrict gws = interface->gws;
  for(int k = 0; k <= ksize+1; k++){
    for(int j = 0; j <= jsize+1; j++){
      for(int i = 1; i <= isize; i++){
        // for (almost) single-phase region,
        //   surface reconstruction is not needed
        if(VOF(i, j, k) < VOFMIN || 1.-VOFMIN < VOF(i, j, k)){
          continue;
        }
        /* ! average nx ! 10 ! */
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
        /* ! average ny ! 10 ! */
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
        /* ! average nz ! 10 ! */
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
        /* ! normalise and obtain center normals ! 12 ! */
        nx = nx / DXF(i);
        ny = ny / dy;
        nz = nz / dz;
        double norm = sqrt(
            +pow(nx, 2.)
            +pow(ny, 2.)
            +pow(nz, 2.)
        );
        norm = fmax(norm, DBL_EPSILON);
        nx /= norm;
        ny /= norm;
        nz /= norm;
        /* Newton-Raphson iteration */
        double a000;
        {
          /* ! compute constants a priori ! 8 ! */
          double exp_psis[ORDER_GAUSS][ORDER_GAUSS][ORDER_GAUSS];
          for(int kk = 0; kk < ORDER_GAUSS; kk++){
            for(int jj = 0; jj < ORDER_GAUSS; jj++){
              for(int ii = 0; ii < ORDER_GAUSS; ii++){
                exp_psis[kk][jj][ii] = exp(-2.*VOFBETA*psi(nx, ny, nz, gps[ii], gps[jj], gps[kk]));
              }
            }
          }
          /* ! initial guess ! 1 ! */
          a000 = 1. / VOF(i, j, k) - 1.;
          int cnt = 0;
          double f, fp;
          do {
            /* ! sum up ! 13 ! */
            f = -VOF(i, j, k);
            fp = 0.;
            for(int kk = 0; kk < ORDER_GAUSS; kk++){
              for(int jj = 0; jj < ORDER_GAUSS; jj++){
                for(int ii = 0; ii < ORDER_GAUSS; ii++){
                  double tmp = gws[ii] * gws[jj] * gws[kk];
                  // function
                  f  += tmp/(1.+exp_psis[kk][jj][ii]*a000);
                  // its derivative (f^{prime})
                  fp -= tmp*exp_psis[kk][jj][ii]/pow(1.+exp_psis[kk][jj][ii]*a000, 2.);
                }
              }
            }
            /* ! update D ! 1 ! */
            a000 -= f/fp;
            cnt++;
            if(cnt > nrcntmax){
              break;
            }
          } while(fabs(f) > nreps);
          /* ! convert D to d ! 1 ! */
          a000 = -1./(2.*VOFBETA)*log(a000);
        }
        /* ! store normal and intercept ! 4 ! */
        NORMAL(i, j, k).a100 = nx;
        NORMAL(i, j, k).a010 = ny;
        NORMAL(i, j, k).a001 = nz;
        NORMAL(i, j, k).a000 = a000;
      }
    }
  }
  return 0;
}

static int compute_curvature(const domain_t *domain, interface_t *interface){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict dxf = domain->dxf;
  const double            dy  = domain->dy;
  const double            dz  = domain->dz;
  const dvof_t * restrict dvof = interface->dvof;
  double * restrict curv = interface->curv;
  for(int k = 0; k <= ksize+1; k++){
    for(int j = 0; j <= jsize+1; j++){
      for(int i = 1; i <= isize; i++){
        /* ! compute mean curvature from corner normals ! 19 ! */
        double dnxdx = 0.25/DXF(i)*(
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

