#include <assert.h>
#include "domain.h"
#include "interface.h"
#include "arrays/interface.h"



int interface_update_boundaries_vof(const domain_t * restrict domain, double * restrict vof){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  /* update y halo values */
  {
    const MPI_Comm comm = domain->sdecomp->comm_cart;
    // check neighbour ranks, negative and positive
    int neg, pos;
    MPI_Cart_shift(comm, 1, 1, &neg, &pos);
    // create datatype
    MPI_Datatype dtype;
    MPI_Type_vector(ksize, 2*(isize+2), (isize+2)*(jsize+4), MPI_DOUBLE, &dtype);
    MPI_Type_commit(&dtype);
    // communicate
    MPI_Sendrecv(
      /* send to   pos. */ &VOF(0, jsize-1, 1), 1, dtype, pos, 0,
      /* recv from neg. */ &VOF(0,      -1, 1), 1, dtype, neg, 0,
      comm, MPI_STATUS_IGNORE
    );
    MPI_Sendrecv(
      /* send to   neg. */ &VOF(0,       1, 1), 1, dtype, neg, 0,
      /* recv from pos. */ &VOF(0, jsize+1, 1), 1, dtype, pos, 0,
      comm, MPI_STATUS_IGNORE
    );
    // clean-up used datatype
    MPI_Type_free(&dtype);
  }
  /* update z halo values */
  {
    const MPI_Comm comm = domain->sdecomp->comm_cart;
    // check neighbour ranks, negative and positive
    int neg, pos;
    MPI_Cart_shift(comm, 2, 1, &neg, &pos);
    // create datatype
    MPI_Datatype dtype;
    MPI_Type_contiguous(
        (isize+2)*(jsize+4)*2,
        MPI_DOUBLE,
        &dtype
    );
    MPI_Type_commit(&dtype);
    // communicate
    MPI_Sendrecv(
      /* send to   pos. */ &VOF(0, -1, ksize-1), 1, dtype, pos, 0,
      /* recv from neg. */ &VOF(0, -1,      -1), 1, dtype, neg, 0,
      comm, MPI_STATUS_IGNORE
    );
    MPI_Sendrecv(
      /* send to   neg. */ &VOF(0, -1,       1), 1, dtype, neg, 0,
      /* recv from pos. */ &VOF(0, -1, ksize+1), 1, dtype, pos, 0,
      comm, MPI_STATUS_IGNORE
    );
    // clean-up used datatype
    MPI_Type_free(&dtype);
  }
  /* ! set boundary values ! 6 ! */
  for(int k = -1; k <= ksize+2; k++){
    for(int j = -1; j <= jsize+2; j++){
      VOF(      0, j, k) = 0.;
      VOF(isize+1, j, k) = 0.;
    }
  }
  return 0;
}

