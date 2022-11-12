#include <string.h>
#include <math.h>
#include "domain.h"
#include "fluid.h"
#include "temperature.h"
#include "interface.h"
#include "arrays/fluid.h"
#include "arrays/temperature.h"
#include "internal.h"


/**
 * @brief comute right-hand-side of Runge-Kutta scheme
 * @param[in   ] domain      : information related to MPI domain decomposition
 * @param[in   ] rkstep      : Runge-Kutta step
 * @param[inout] fluid       : velocity and pressure (in), RK source terms (inout)
 * @param[in   ] temperature : buoyancy force (in)
 * @return                   : error code
 */
int fluid_compute_rhs(const domain_t * restrict domain, const int rkstep, fluid_t * restrict fluid, const temperature_t * restrict temperature, const interface_t * restrict interface){
  compute_rhs_ux(domain, rkstep, fluid, temperature, interface);
  compute_rhs_uy(domain, rkstep, fluid, interface);
#if NDIMS == 3
  compute_rhs_uz(domain, rkstep, fluid, interface);
#endif
  return 0;
}

