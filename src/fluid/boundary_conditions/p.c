#include "domain.h"
#include "fluid.h"
#include "arrays/fluid.h"



/**
 * @brief update boundary values of pressure and scalar potential
 * @param[in   ] domain : information about domain decomposition and size
 * @param[inout] p      : pressure or scalar potential
 * @return              : error code
 */
int fluid_update_boundaries_p(const domain_t * restrict domain, double * restrict p){
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
      /* send to   pos. */ &P(1,   jsize, 1), 1, dtype, pos, 0,
      /* recv from neg. */ &P(1,       0, 1), 1, dtype, neg, 0,
      comm, MPI_STATUS_IGNORE
    );
    MPI_Sendrecv(
      /* send to   neg. */ &P(1,       1, 1), 1, dtype, neg, 0,
      /* recv from pos. */ &P(1, jsize+1, 1), 1, dtype, pos, 0,
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
      /* send to   pos. */ &P(1, 1,   ksize), 1, dtype, pos, 0,
      /* recv from neg. */ &P(1, 1,       0), 1, dtype, neg, 0,
      comm, MPI_STATUS_IGNORE
    );
    MPI_Sendrecv(
      /* send to   neg. */ &P(1, 1,       1), 1, dtype, neg, 0,
      /* recv from pos. */ &P(1, 1, ksize+1), 1, dtype, pos, 0,
      comm, MPI_STATUS_IGNORE
    );
    // clean-up used datatype
    MPI_Type_free(&dtype);
  }
  /* ! set boundary values ! 6 ! */
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      P(      0, j, k) = P(    1, j, k); // Neumann
      P(isize+1, j, k) = P(isize, j, k); // Neumann
    }
  }
  return 0;
}

