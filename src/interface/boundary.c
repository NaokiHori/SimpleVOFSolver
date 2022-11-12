#include <assert.h>
#include "domain.h"
#include "interface.h"
#include "arrays/interface.h"


#if NDIMS == 2

int interface_update_boundaries_vof(const domain_t *domain, double *vof){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  {
    const MPI_Comm comm = domain->sdecomp->comm_cart;
    // check neighbour ranks, negative and positive
    int ymrank, yprank;
    MPI_Cart_shift(comm, 1, 1, &ymrank, &yprank);
    // create datatype
    MPI_Datatype dtype;
    MPI_Type_contiguous(2*(isize+2), MPI_DOUBLE, &dtype);
    MPI_Type_commit(&dtype);
    // send in positive direction
    MPI_Sendrecv(
      &VOF(0,   jsize-1), 1, dtype, yprank, 0,
      &VOF(0,        -1), 1, dtype, ymrank, 0,
      comm, MPI_STATUS_IGNORE
    );
    // send in negative direction
    MPI_Sendrecv(
      &VOF(0,         1), 1, dtype, ymrank, 0,
      &VOF(0,   jsize+1), 1, dtype, yprank, 0,
      comm, MPI_STATUS_IGNORE
    );
    // clean-up used datatype
    MPI_Type_free(&dtype);
  }
  for(int j = -1; j <= jsize+2; j++){
    VOF(      0, j) = 0.;
    VOF(isize+1, j) = 0.;
  }
  return 0;
}

#else // NDIMS == 3

int interface_update_boundaries_vof(const domain_t *domain, double *vof){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  /* update y halo values */
  {
    const MPI_Comm comm = domain->sdecomp->comm_cart;
    // check neighbour ranks, negative and positive
    int ymrank, yprank;
    MPI_Cart_shift(comm, 1, 1, &ymrank, &yprank);
    // create datatype
    MPI_Datatype dtype;
    MPI_Type_vector(
        ksize+4,
        2*(isize+2),
        (isize+2)*(jsize+4),
        MPI_DOUBLE,
        &dtype
    );
    MPI_Type_commit(&dtype);
    // send in positive direction
    MPI_Sendrecv(
      &VOF(0,   jsize-1, -1), 1, dtype, yprank, 0,
      &VOF(0,        -1, -1), 1, dtype, ymrank, 0,
      comm, MPI_STATUS_IGNORE
    );
    // send in negative direction
    MPI_Sendrecv(
      &VOF(0,         1, -1), 1, dtype, ymrank, 0,
      &VOF(0,   jsize+1, -1), 1, dtype, yprank, 0,
      comm, MPI_STATUS_IGNORE
    );
    // clean-up used datatype
    MPI_Type_free(&dtype);
  }
  /* update z halo values */
  {
    const MPI_Comm comm = domain->sdecomp->comm_cart;
    // check neighbour ranks, negative and positive
    int zmrank, zprank;
    MPI_Cart_shift(comm, 2, 1, &zmrank, &zprank);
    // create datatype
    MPI_Datatype dtype;
    MPI_Type_contiguous(
        (isize+2)*(jsize+4)*2,
        MPI_DOUBLE,
        &dtype
    );
    MPI_Type_commit(&dtype);
    // send in positive direction
    MPI_Sendrecv(
      &VOF(0, -1, ksize-1), 1, dtype, zprank, 0,
      &VOF(0, -1,      -1), 1, dtype, zmrank, 0,
      comm, MPI_STATUS_IGNORE
    );
    // send in negative direction
    MPI_Sendrecv(
      &VOF(0, -1,       1), 1, dtype, zmrank, 0,
      &VOF(0, -1, ksize+1), 1, dtype, zprank, 0,
      comm, MPI_STATUS_IGNORE
    );
    // clean-up used datatype
    MPI_Type_free(&dtype);
  }
  for(int k = -1; k <= ksize+2; k++){
    for(int j = -1; j <= jsize+2; j++){
      VOF(      0, j, k) = 0.;
      VOF(isize+1, j, k) = 0.;
    }
  }
  return 0;
}

#endif // NDIMS
