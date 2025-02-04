#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid_solver.h"
#include "array_macros/fluid/psi.h"

static int assign_boundary_conditions_in_x(
    const domain_t * domain,
    double * psi
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  // set boundary values
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      PSI(      0, j, k) = PSI(    1, j, k); // Neumann
      PSI(isize+1, j, k) = PSI(isize, j, k); // Neumann
    }
  }
  return 0;
}

/**
 * @brief update boundary values of the scalar potential
 * @param[in]     domain : information about domain decomposition and size
 * @param[in,out] p      : pressure
 * @return               : error code
 */
int fluid_update_boundaries_psi(
    const domain_t * domain,
    array_t * psi
){
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
    MPI_DOUBLE,
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, psi)){
    return 1;
  }
  if(0 != halo_communicate_in_z(domain, dtypes + 1, psi)){
    return 1;
  }
  assign_boundary_conditions_in_x(domain, psi->data);
  return 0;
}

