#include <math.h>
#include "interface.h"
#include "internal.h"


#if NDIMS == 2

double psi(const double a10, const double a01, const double x, const double y){
  return
    +a10*x
    +a01*y
  ;
}

double H(const double a00, const double a10, const double a01, const double x, const double y){
  return 1./(1.+exp(-2.*VOFBETA*(psi(a10, a01, x, y)+a00)));
}

#else // NDIMS == 3

double psi(const double a100, const double a010, const double a001, const double x, const double y, const double z){
  return
    +a100*x
    +a010*y
    +a001*z
  ;
}

double H(const double a000, const double a100, const double a010, const double a001, const double x, const double y, const double z){
  return 1./(1.+exp(-2.*VOFBETA*(psi(a100, a010, a001, x, y, z)+a000)));
}

#endif // NDIMS
