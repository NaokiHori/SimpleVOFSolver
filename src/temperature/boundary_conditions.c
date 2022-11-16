#include "domain.h"
#include "temperature.h"
#include "arrays/temperature.h"


#if NDIMS == 2

/**
 * @brief update boundary values of temperature
 * @param[in   ] domain : information about domain decomposition and size
 * @param[inout] temp   : temperature
 * @return              : error code
 */
int temperature_update_boundaries_temp(const domain_t * restrict domain, double * restrict temp){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  /* ! update y halo values in 2D ! 23 ! */
  {
    const MPI_Comm comm = domain->sdecomp->comm_cart;
    // check neighbour ranks, negative and positive
    int neg, pos;
    MPI_Cart_shift(comm, 1, 1, &neg, &pos);
    // create datatype
    MPI_Datatype dtype;
    MPI_Type_contiguous(isize, MPI_DOUBLE, &dtype);
    MPI_Type_commit(&dtype);
    // communicate
    MPI_Sendrecv(
      /* send to   pos. */ &TEMP(1,   jsize), 1, dtype, pos, 0,
      /* recv from neg. */ &TEMP(1,       0), 1, dtype, neg, 0,
      comm, MPI_STATUS_IGNORE
    );
    MPI_Sendrecv(
      /* send to   neg. */ &TEMP(1,       1), 1, dtype, neg, 0,
      /* recv from pos. */ &TEMP(1, jsize+1), 1, dtype, pos, 0,
      comm, MPI_STATUS_IGNORE
    );
    // clean-up used datatype
    MPI_Type_free(&dtype);
  }
  /* ! set boundary values ! 10 ! */
  // values do not matter,
  //   e.g., combination of temp_xm = 1, temp_xp = 0 also works,
  //   but "temp_xm - temp_xp" should be 1,
  //   since governing equations are normalised with this assumption
  const double temp_xm = +0.5;
  const double temp_xp = -0.5;
  for(int j = 1; j <= jsize; j++){
    TEMP(      0, j) = temp_xm;
    TEMP(isize+1, j) = temp_xp;
  }
  return 0;
}

#else // NDIMS == 3

/**
 * @brief update boundary values of temperature
 * @param[in   ] domain : information about domain decomposition and size
 * @param[inout] temp   : temperature
 * @return              : error code
 */
int temperature_update_boundaries_temp(const domain_t * restrict domain, double * restrict temp){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  /* ! update y halo values in 3D ! 23 ! */
  {
    const MPI_Comm comm = domain->sdecomp->comm_cart;
    // check neighbour ranks, negative and positive
    int neg, pos;
    MPI_Cart_shift(comm, 1, 1, &neg, &pos);
    // create datatype
    MPI_Datatype dtype;
    MPI_Type_vector(ksize, isize, (isize+2)*(jsize+2), MPI_DOUBLE, &dtype);
    MPI_Type_commit(&dtype);
    // communicate
    MPI_Sendrecv(
      /* send to   pos. */ &TEMP(1,   jsize, 1), 1, dtype, pos, 0,
      /* recv from neg. */ &TEMP(1,       0, 1), 1, dtype, neg, 0,
      comm, MPI_STATUS_IGNORE
    );
    MPI_Sendrecv(
      /* send to   neg. */ &TEMP(1,       1, 1), 1, dtype, neg, 0,
      /* recv from pos. */ &TEMP(1, jsize+1, 1), 1, dtype, pos, 0,
      comm, MPI_STATUS_IGNORE
    );
    // clean-up used datatype
    MPI_Type_free(&dtype);
  }
  /* ! update z halo values in 3D ! 23 ! */
  {
    const MPI_Comm comm = domain->sdecomp->comm_cart;
    // check neighbour ranks, negative and positive
    int neg, pos;
    MPI_Cart_shift(comm, 2, 1, &neg, &pos);
    // create datatype
    MPI_Datatype dtype;
    MPI_Type_vector(jsize, isize, isize+2, MPI_DOUBLE, &dtype);
    MPI_Type_commit(&dtype);
    // communicate
    MPI_Sendrecv(
      /* send to   pos. */ &TEMP(1, 1,   ksize), 1, dtype, pos, 0,
      /* recv from neg. */ &TEMP(1, 1,       0), 1, dtype, neg, 0,
      comm, MPI_STATUS_IGNORE
    );
    MPI_Sendrecv(
      /* send to   neg. */ &TEMP(1, 1,       1), 1, dtype, neg, 0,
      /* recv from pos. */ &TEMP(1, 1, ksize+1), 1, dtype, pos, 0,
      comm, MPI_STATUS_IGNORE
    );
    // clean-up used datatype
    MPI_Type_free(&dtype);
  }
  /* ! set boundary values ! 12 ! */
  // values do not matter,
  //   e.g., combination of temp_xm = 1, temp_xp = 0 also works,
  //   but "temp_xm - temp_xp" should be 1,
  //   since governing equations are normalised with this assumption
  const double temp_xm = +0.5;
  const double temp_xp = -0.5;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      TEMP(      0, j, k) = temp_xm;
      TEMP(isize+1, j, k) = temp_xp;
    }
  }
  return 0;
}

#endif // NDIMS
