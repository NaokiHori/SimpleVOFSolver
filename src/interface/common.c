#include <math.h>
#include "interface.h"
#include "internal.h"



double psi(const double a10, const double a01, const double x, const double y){
  return
    +a10*x
    +a01*y
  ;
}

double H(const double a00, const double a10, const double a01, const double x, const double y){
  return 1./(1.+exp(-2.*VOFBETA*(psi(a10, a01, x, y)+a00)));
}

