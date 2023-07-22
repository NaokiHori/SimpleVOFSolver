#include "array.h"
#include "domain.h"
#include "halo.h"
#include "interface_solver.h"
#include "array_macros/interface/vof.h"

static int assign_boundary_conditions_in_x(
    const domain_t * domain,
    double * vof
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  // set boundary values
#if NDIMS == 2
  for(int j = -1; j <= jsize + 2; j++){
    VOF(      0, j) = 0.;
    VOF(isize+1, j) = 0.;
  }
#else
  for(int k = -1; k <= ksize + 2; k++){
    for(int j = -1; j <= jsize + 2; j++){
      VOF(      0, j, k) = 0.;
      VOF(isize+1, j, k) = 0.;
    }
  }
#endif
  return 0;
}

/**
 * @brief update boundary values of vof field
 * @param[in]     domain : information about domain decomposition and size
 * @param[in,out] vof    : indicator function
 * @return               : error code
 */
int interface_update_boundaries_vof(
    const domain_t * domain,
    array_t * vof
){
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
#if NDIMS == 3
    MPI_DOUBLE,
#endif
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, vof)){
    return 1;
  }
#if NDIMS == 3
  if(0 != halo_communicate_in_z(domain, dtypes + 1, vof)){
    return 1;
  }
#endif
  assign_boundary_conditions_in_x(domain, vof->data);
  return 0;
}

