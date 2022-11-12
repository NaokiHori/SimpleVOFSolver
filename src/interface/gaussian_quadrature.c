#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "common.h"


static int get_f(const int n, double *f0){
  /* Returns coefficients of N-th order Legendre polynomial */
  /* by using the following recursive relation */
  /* L(N, x) = (2n-1)/n*L(N-1, x) - (n-1)/n*L(N-2, x) */
  /* REFERENCE: https://en.wikipedia.org/wiki/Legendre_polynomials#Recurrence_relations */
  if(n < 0){
    printf("ORDER_GAUSS is negative: %3d\n", n);
    exit(EXIT_FAILURE);
    return 1;
  }else if(n == 0){
    f0[0] = 1.;
  }else if(n == 1){
    f0[0] = 0.;
    f0[1] = 1.;
  }else{
    double *f1 = NULL;
    double *f2 = NULL;
    f1 = common_calloc(n  , sizeof(double));
    f2 = common_calloc(n-1, sizeof(double));
    get_f(n-1, f1);
    get_f(n-2, f2);
    /* L(N, x) = (2n-1)/n*L(N-1, x) - (n-1)/n*L(N-2, x) */
    // first multiply pre-factors
    for(int i=0; i<n; i++){
      f1[i] *= (2.*n-1.)/(1.*n);
    }
    for(int i=0; i<n-1; i++){
      f2[i] *= -(1.*n-1.)/(1.*n);
    }
    // sum up N-1 and N-2 -th order polynomials (f1 and f2, respectively)
    for(int i=0; i<n+1; i++){
      if(i == 0){
        f0[i] = 0.+f2[i];
      }else if(i < n-1){
        f0[i] = f1[i-1]+f2[i];
      }else{
        f0[i] = f1[i-1]+0.;
      }
    }
    common_free(f1);
    common_free(f2);
  }
  return 0;
}

static int get_fp(const int n, const double *f, double *fp){
  /* compute derivative of f: f^{\prime} */
  /* we readily find the relation */
  /* (a_n x^n)^{\prime} = a_n n x^{n-1} */
  for(int i=1; i<n+1; i++){
    fp[i-1] = 1.*i*f[i];
  }
  return 0;
}

static double compute_polynomial(const int n, const double *f, const double x){
  /* return f(N, x): N-th order polynomial at position x whose coefficients are f */
  double val = 0.;
  for(int i=0; i<n+1; i++){
    val += f[i]*pow(x, 1.*i);
  }
  return val;
}

static double find_roots(const int n, const double *f, const double *fp, double *roots, const double epsilon, const int maxiter){
  /* find roots of the given Legendre polynomial f */
  /* return the maximum residual to evaluate the quality */
  // set ICs
  // in order to avoid the same answer
  // each IC should be close enough to the answer
  const double dx = (1.-(-1.))/pow(100*n, 2.);
  double lx, rx, lval, rval;
  double max_residual;
  lx = -1.;
  lval = compute_polynomial(n, f, lx);
  for(int i=0; i<n; i++){
    for(;;){
      rx = lx+dx;
      rval = compute_polynomial(n, f, rx);
      if(lval*rval < 0.){
        roots[i] = 0.5*(lx+rx);
        lx = rx;
        lval = rval;
        break;
      }else{
        lx = rx;
        lval = rval;
      }
    }
  }
  // main part to find roots
  max_residual = 0.;
  for(int i=0; i<n; i++){
    int cnt;
    double residual;
    cnt = 0;
    while(++cnt){
      double root = roots[i];
      double local_f  = compute_polynomial(n,   f,  root);
      double local_fp = compute_polynomial(n-1, fp, root);
      double mul;
      mul = 1.;
      for(int j=0; j<i; j++){
        mul *= 1./(root-roots[j]);
      }
      root -= local_f/(local_fp-local_f*mul);
      residual = fabs(root-roots[i]);
      if(residual < epsilon){
        max_residual = fmax(residual, max_residual);
        break;
      }
      if(cnt > maxiter){
        max_residual = fmax(residual, max_residual);
        break;
      }
      roots[i] = root;
    }
  }
  return max_residual;
}

static int find_weights(const int n, const double *fp, const double *roots, double *weights){
  for(int i=0; i<n; i++){
    double x = roots[i];
    weights[i] = 2./(1.-pow(x, 2.))/pow(compute_polynomial(n-1, fp, x), 2.);
  }
  return 0;
}


static double generate_gaussian_quadrature(const int n, double *roots, double *weights){
  double *f = NULL;
  double *fp = NULL;
  double max_residual;
  // N-th order polynomial has N+1 coefficients
  f  = common_calloc(n+1, sizeof(double));
  // its derivative has N coefficients
  fp = common_calloc(n  , sizeof(double));
  // first compute coefficients of N-th order polynomial and its derivative
  get_f(n, f);
  get_fp(n, f, fp);
  // find roots
  max_residual = find_roots(n, f, fp, roots, DBL_EPSILON, 10000);
  // and weights
  find_weights(n, fp, roots, weights);
  //
  common_free(f);
  common_free(fp);
  return max_residual;
}

double init_gaussian_quadrature(const int n, double *gps, double *gws){
  /*
   * computes Gaussian quadrature (roots: gps, weights: gws) with interval [-1, 1]
   * returns maximum residual, which should be small enough
   * in order to change interval, use "convert_gaussian_quadrature"
   */
  if(n > 6){
    printf("WARNING: ORDER_GAUSS is too large: %3d\n", n);
    printf("This would not improve the accuracy of VOF algorithm\n");
    printf("Also the generated Gaussian quadrature can be inaccurate\n");
  }
  // compute roots (gps) and weights (gws), return residual
  return generate_gaussian_quadrature(n, gps, gws);
}

int convert_gaussian_quadrature(const int n, const double xm, const double xp, const double *gps_bef, const double *gws_bef, double *gps_aft, double *gws_aft){
  /*
   * change interval from [-1, 1] to the given one [xm, xp]
   * see https://en.wikipedia.org/wiki/Gaussian_quadrature#Change_of_interval
   */
  for(int i=0; i<n; i++){
    gps_aft[i] = 0.5*(xp-xm)*gps_bef[i]+0.5*(xp+xm);
    gws_aft[i] = 0.5*(xp-xm)*gws_bef[i];
  }
  return 0;
}

