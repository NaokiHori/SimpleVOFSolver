#include <stdbool.h>
#include "runge_kutta.h"
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "interface.h"
#include "interface_solver.h"
#include "integrate.h"

// integrate the equations for one time step
int integrate(
    const domain_t * domain,
    fluid_t * fluid,
    interface_t * interface,
    double * dt
){
  // check grid in x direction is uniform
  //   if it is, I can use more efficient algorithm
  //   to compute the scalar potential;
  //   otherwise a versatile version is adopted
  bool x_grid_is_uniform = false;
  domain_check_x_grid_is_uniform(domain, &x_grid_is_uniform);
  // decide time step size
  if(0 != fluid_decide_dt(domain, fluid, dt)){
    return 1;
  }
  // Runge-Kutta iterations
  // max iteration, should be three
  const size_t rkstepmax = sizeof(rkcoefs) / sizeof(rkcoef_t);
  for(size_t rkstep = 0; rkstep < rkstepmax; rkstep++){
    // update vof field
    if(0 != interface_compute_curvature_tensor(domain, interface)){
      return 1;
    }
    if(0 != interface_compute_force(domain, interface)){
      return 1;
    }
    if(0 != interface_update_vof(domain, rkstep, *dt, fluid, interface)){
      return 1;
    }
    // predict flow field
    if(0 != fluid_predict_field(domain, rkstep, *dt, fluid, interface)){
      return 1;
    }
    // now the temperature field has been updated,
    //   while the velocity field is not divergence free
    //   and thus the following correction step is needed
    // compute scalar potential
    if(x_grid_is_uniform){
      if(0 != fluid_compute_potential_dct(domain, rkstep, *dt, fluid)){
        return 1;
      }
    }else{
      if(0 != fluid_compute_potential_dft(domain, rkstep, *dt, fluid)){
        return 1;
      }
    }
    // correct velocity field to satisfy mass conservation
    if(0 != fluid_correct_velocity(domain, rkstep, *dt, fluid)){
      return 1;
    }
    // update pressure
    if(0 != fluid_update_pressure(domain, rkstep, *dt, fluid)){
      return 1;
    }
  }
  return 0;
}

