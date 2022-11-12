#include <stdio.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include "domain.h"
#include "fileio.h"
#include "interface.h"
#include "arrays/interface.h"
#include "internal.h"


#if NDIMS == 2

int check_interface(const char fname[], const domain_t *domain, const double time, const interface_t *interface){
  const MPI_Comm comm_cart = domain->sdecomp->comm_cart;
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double dx = domain->dx;
  const double dy = domain->dy;
  const double *vof = interface->vof;
  double min = DBL_MAX;
  double max = DBL_MIN;
  double sum = 0.;
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      min = fmin(min, VOF(i, j));
      max = fmax(max, VOF(i, j));
      sum += VOF(i, j)*dx*dy;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &min, 1, MPI_DOUBLE, MPI_MIN, comm_cart);
  MPI_Allreduce(MPI_IN_PLACE, &max, 1, MPI_DOUBLE, MPI_MAX, comm_cart);
  MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, comm_cart);
  int myrank;
  MPI_Comm_rank(domain->sdecomp->comm_cart, &myrank);
  if(myrank == 0){
    FILE *fp = fileio_fopen(fname, "a");
    if(fp != NULL){
      fprintf(fp, "%8.2f % 18.15e % 18.15e % 18.15e\n", time, min, max, sum);
      fileio_fclose(fp);
    }
  }
  return 0;
}

#else // NDIMS == 3

int check_interface(const char fname[], const domain_t *domain, const double time, const interface_t *interface){
  const MPI_Comm comm_cart = domain->sdecomp->comm_cart;
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double dx = domain->dx;
  const double dy = domain->dy;
  const double dz = domain->dz;
  const double *vof = interface->vof;
  double min = DBL_MAX;
  double max = DBL_MIN;
  double sum = 0.;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        min = fmin(min, VOF(i, j, k));
        max = fmax(max, VOF(i, j, k));
        sum += VOF(i, j, k)*dx*dy*dz;
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &min, 1, MPI_DOUBLE, MPI_MIN, comm_cart);
  MPI_Allreduce(MPI_IN_PLACE, &max, 1, MPI_DOUBLE, MPI_MAX, comm_cart);
  MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, comm_cart);
  int myrank;
  MPI_Comm_rank(domain->sdecomp->comm_cart, &myrank);
  if(myrank == 0){
    FILE *fp = fileio_fopen(fname, "a");
    if(fp != NULL){
      fprintf(fp, "%8.2f % 18.15e % 18.15e % 18.15e\n", time, min, max, sum);
      fileio_fclose(fp);
    }
  }
  return 0;
}

#endif // NDIMS
