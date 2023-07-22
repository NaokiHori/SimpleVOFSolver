#include <stdio.h>
#include <math.h>
#include "domain.h"
#include "interface.h"
#include "fileio.h"
#include "array_macros/domain/dxf.h"
#include "array_macros/interface/vof.h"
#include "internal.h"

int logging_check_vof(
    const char fname[],
    const domain_t * domain,
    const double time,
    const interface_t * interface
){
  const int root = 0;
  int myrank = root;
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_rank(domain->info, &myrank);
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * dxf = domain->dxf;
  const double dy = domain->dy;
#if NDIMS == 3
  const double dz = domain->dz;
#endif
  const double * vof = interface->vof.data;
  double min = 1.;
  double max = 0.;
  double sums[2] = {0.};
#if NDIMS == 2
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      const double dx = DXF(i  );
      const double cellsize = dx * dy;
      const double lvof = VOF(i, j);
      min = fmin(min, lvof);
      max = fmax(max, lvof);
      sums[0] += lvof * cellsize;
      sums[1] +=        cellsize;
    }
  }
#else
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        const double dx = DXF(i  );
        const double cellsize = dx * dy * dz;
        const double lvof = VOF(i, j, k);
        min = fmin(min, lvof);
        max = fmax(max, lvof);
        sums[0] += lvof * cellsize;
        sums[1] +=        cellsize;
      }
    }
  }
#endif
  MPI_Allreduce(MPI_IN_PLACE, &min, 1, MPI_DOUBLE, MPI_MIN, comm_cart);
  MPI_Allreduce(MPI_IN_PLACE, &max, 1, MPI_DOUBLE, MPI_MAX, comm_cart);
  MPI_Allreduce(MPI_IN_PLACE, sums, 2, MPI_DOUBLE, MPI_SUM, comm_cart);
  if(root == myrank){
    FILE * fp = fileio.fopen(fname, "a");
    if(NULL == fp){
      return 0;
    }
    fprintf(fp, "%8.2f ", time);
    fprintf(fp, "% 18.15e ", fabs(min) < 1.e-100 ? 1.e-99 : min);
    fprintf(fp, "% 18.15e ", max);
    fprintf(fp, "% 18.15e\n", sums[0] / sums[1]);
    fileio.fclose(fp);
  }
  return 0;
}

