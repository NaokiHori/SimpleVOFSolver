#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "common.h"
#include "config.h"
#include "domain.h"
#include "interface.h"
#include "fileio.h"
#include "internal.h"
#include "arrays/interface.h"


#if NDIMS == 2

static interface_t *allocate(const domain_t *domain){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  interface_t *interface = common_calloc(1, sizeof(interface_t));
  interface->vof      = common_calloc(VOF_SIZE_0      * VOF_SIZE_1     , sizeof(double));
  interface->dvof     = common_calloc(DVOF_SIZE_0     * DVOF_SIZE_1    , sizeof(dvof_t));
  interface->normal   = common_calloc(NORMAL_SIZE_0   * NORMAL_SIZE_1  , sizeof(nrml_t));
  interface->curv     = common_calloc(CURV_SIZE_0     * CURV_SIZE_1    , sizeof(double));
  interface->voffluxx = common_calloc(VOFFLUXX_SIZE_0 * VOFFLUXX_SIZE_1, sizeof(double));
  interface->voffluxy = common_calloc(VOFFLUXY_SIZE_0 * VOFFLUXY_SIZE_1, sizeof(double));
  interface->voffrcx  = common_calloc(VOFFRCX_SIZE_0  * VOFFRCX_SIZE_1 , sizeof(double));
  interface->voffrcy  = common_calloc(VOFFRCY_SIZE_0  * VOFFRCY_SIZE_1 , sizeof(double));
  interface->vofsrca  = common_calloc(VOFSRCA_SIZE_0  * VOFSRCA_SIZE_1 , sizeof(double));
  interface->vofsrcb  = common_calloc(VOFSRCB_SIZE_0  * VOFSRCB_SIZE_1 , sizeof(double));
  interface->gps      = common_calloc(ORDER_GAUSS, sizeof(double));
  interface->gws      = common_calloc(ORDER_GAUSS, sizeof(double));
  return interface;
}

static int init(const domain_t *domain, interface_t *interface){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ioffs = domain->offsets[0];
  const int joffs = domain->offsets[1];
  const double lx = domain->lengths[0];
  const double ly = domain->lengths[1];
  const double dx = domain->dx;
  const double dy = domain->dy;
  const double cx = 0.50*lx;
  const double cy = 0.50*ly;
  const double delta = fmin(dx, dy);
  double xoffs = dx * ioffs;
  double yoffs = dy * joffs;
  double *vof = interface->vof;
  double *g_gps = NULL;
  double *g_gws = NULL;
  double *l_gps_x = NULL;
  double *l_gws_x = NULL;
  double *l_gps_y = NULL;
  double *l_gws_y = NULL;
  g_gps   = common_calloc(ORDER_GAUSS, sizeof(double));
  g_gws   = common_calloc(ORDER_GAUSS, sizeof(double));
  l_gps_x = common_calloc(ORDER_GAUSS, sizeof(double));
  l_gws_x = common_calloc(ORDER_GAUSS, sizeof(double));
  l_gps_y = common_calloc(ORDER_GAUSS, sizeof(double));
  l_gws_y = common_calloc(ORDER_GAUSS, sizeof(double));
  // obtain gaussian quadrature variables with interval [-1, 1]
  // (interface->g[pw]s cannot be used since they are with [-0.5, 0.5])
  init_gaussian_quadrature(ORDER_GAUSS, g_gps, g_gws);
  for(int j = 1; j <= jsize; j++){
    double y0 = 0.5*(2*j-1)*dy+yoffs;
    // compute local gaussian quadrature with interval [y-0.5*dy, y+0.5*dy]
    convert_gaussian_quadrature(ORDER_GAUSS, y0-0.5*dy, y0+0.5*dy, g_gps, g_gws, l_gps_y, l_gws_y);
    for(int i = 1; i <= isize; i++){
      double x0 = 0.5*(2*i-1)*dx+xoffs;
      // compute local gaussian quadrature with interval [x-0.5*dx, x+0.5*dx]
      convert_gaussian_quadrature(ORDER_GAUSS, x0-0.5*dx, x0+0.5*dx, g_gps, g_gws, l_gps_x, l_gws_x);
      // integrate H(distance) to determine VOF value
      for(int jj=0; jj<ORDER_GAUSS; jj++){
        for(int ii=0; ii<ORDER_GAUSS; ii++){
          double dist = l_gps_y[jj]-cy;
          VOF(i, j) +=
            l_gws_x[ii]*l_gws_y[jj]
            *0.5*(1.+tanh(VOFBETA*dist/delta))
            /dx/dy;
        }
      }
    }
  }
  common_free(g_gps);
  common_free(g_gws);
  common_free(l_gps_x);
  common_free(l_gws_x);
  common_free(l_gps_y);
  common_free(l_gws_y);
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      VOF(i, j) = fmin(1., VOF(i, j));
      VOF(i, j) = fmax(0., VOF(i, j));
    }
  }
  interface_update_boundaries_vof(domain, vof);
  return 0;
}

static int load(const domain_t *domain, interface_t *interface){
  const int glisize = domain->glsizes[0];
  const int gljsize = domain->glsizes[1];
  const int   isize = domain->mysizes[0];
  const int   jsize = domain->mysizes[1];
  const int ioffset = domain->offsets[0];
  const int joffset = domain->offsets[1];
  // vof: [0:isize+1] x [1:jsize] x [1:ksize]
  double * restrict vof = interface->vof;
  const int glsizes[NDIMS] = {gljsize, glisize+2};
  const int mysizes[NDIMS] = {  jsize,   isize+2};
  const int offsets[NDIMS] = {joffset, ioffset  };
  double *buf = common_calloc(mysizes[0] * mysizes[1], sizeof(double));
  const char *dirname = config.get_string("restart_dir");
  fileio_r_nd_parallel(dirname, "vof", NDIMS, glsizes, mysizes, offsets, buf);
  for(int cnt = 0, j = 1; j <= jsize; j++){
    for(int i = 0; i <= isize+1; i++){
      VOF(i, j) = buf[cnt];
      cnt++;
    }
  }
  common_free(buf);
  interface_update_boundaries_vof(domain, vof);
  return 0;
}

#else // NDIMS == 3

static interface_t *allocate(const domain_t *domain){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  interface_t *interface = common_calloc(1, sizeof(interface_t));
  interface->vof      = common_calloc(VOF_SIZE_0      * VOF_SIZE_1     * VOF_SIZE_2     , sizeof(double));
  interface->dvof     = common_calloc(DVOF_SIZE_0     * DVOF_SIZE_1    * DVOF_SIZE_2    , sizeof(dvof_t));
  interface->normal   = common_calloc(NORMAL_SIZE_0   * NORMAL_SIZE_1  * NORMAL_SIZE_2  , sizeof(nrml_t));
  interface->curv     = common_calloc(CURV_SIZE_0     * CURV_SIZE_1    * CURV_SIZE_2    , sizeof(double));
  interface->voffluxx = common_calloc(VOFFLUXX_SIZE_0 * VOFFLUXX_SIZE_1* VOFFLUXX_SIZE_2, sizeof(double));
  interface->voffluxy = common_calloc(VOFFLUXY_SIZE_0 * VOFFLUXY_SIZE_1* VOFFLUXY_SIZE_2, sizeof(double));
  interface->voffluxz = common_calloc(VOFFLUXZ_SIZE_0 * VOFFLUXZ_SIZE_1* VOFFLUXZ_SIZE_2, sizeof(double));
  interface->voffrcx  = common_calloc(VOFFRCX_SIZE_0  * VOFFRCX_SIZE_1 * VOFFRCX_SIZE_2 , sizeof(double));
  interface->voffrcy  = common_calloc(VOFFRCY_SIZE_0  * VOFFRCY_SIZE_1 * VOFFRCY_SIZE_2 , sizeof(double));
  interface->voffrcz  = common_calloc(VOFFRCZ_SIZE_0  * VOFFRCZ_SIZE_1 * VOFFRCZ_SIZE_2 , sizeof(double));
  interface->vofsrca  = common_calloc(VOFSRCA_SIZE_0  * VOFSRCA_SIZE_1 * VOFSRCA_SIZE_2 , sizeof(double));
  interface->vofsrcb  = common_calloc(VOFSRCB_SIZE_0  * VOFSRCB_SIZE_1 * VOFSRCB_SIZE_2 , sizeof(double));
  interface->gps      = common_calloc(ORDER_GAUSS, sizeof(double));
  interface->gws      = common_calloc(ORDER_GAUSS, sizeof(double));
  return interface;
}

static int init(const domain_t *domain, interface_t *interface){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const int ioffs = domain->offsets[0];
  const int joffs = domain->offsets[1];
  const int koffs = domain->offsets[2];
  const double ly = domain->lengths[1];
  const double dx = domain->dx;
  const double dy = domain->dy;
  const double dz = domain->dz;
  const double delta = fmin(dx, fmin(dy, dz));
  double xoffs = dx * ioffs;
  double yoffs = dy * joffs;
  double zoffs = dz * koffs;
  double *vof = interface->vof;
  double *g_gps = NULL;
  double *g_gws = NULL;
  double *l_gps_x = NULL;
  double *l_gws_x = NULL;
  double *l_gps_y = NULL;
  double *l_gws_y = NULL;
  double *l_gps_z = NULL;
  double *l_gws_z = NULL;
  g_gps   = common_calloc(ORDER_GAUSS, sizeof(double));
  g_gws   = common_calloc(ORDER_GAUSS, sizeof(double));
  l_gps_x = common_calloc(ORDER_GAUSS, sizeof(double));
  l_gws_x = common_calloc(ORDER_GAUSS, sizeof(double));
  l_gps_y = common_calloc(ORDER_GAUSS, sizeof(double));
  l_gws_y = common_calloc(ORDER_GAUSS, sizeof(double));
  l_gps_z = common_calloc(ORDER_GAUSS, sizeof(double));
  l_gws_z = common_calloc(ORDER_GAUSS, sizeof(double));
  const double cy = 0.50*ly;
  // obtain gaussian quadrature variables with interval [-1, 1]
  // (interface->g[pw]s cannot be used since they are with [-0.5, 0.5])
  init_gaussian_quadrature(ORDER_GAUSS, g_gps, g_gws);
  for(int k = 1; k <= ksize; k++){
    double z0 = 0.5*(2*k-1)*dz+zoffs;
    // compute local gaussian quadrature with interval [z-0.5*dz, z+0.5*dz]
    convert_gaussian_quadrature(ORDER_GAUSS, z0-0.5*dz, z0+0.5*dz, g_gps, g_gws, l_gps_z, l_gws_z);
    for(int j = 1; j <= jsize; j++){
      double y0 = 0.5*(2*j-1)*dy+yoffs;
      // compute local gaussian quadrature with interval [y-0.5*dy, y+0.5*dy]
      convert_gaussian_quadrature(ORDER_GAUSS, y0-0.5*dy, y0+0.5*dy, g_gps, g_gws, l_gps_y, l_gws_y);
      for(int i = 1; i <= isize; i++){
        double x0 = 0.5*(2*i-1)*dx+xoffs;
        // compute local gaussian quadrature with interval [x-0.5*dx, x+0.5*dx]
        convert_gaussian_quadrature(ORDER_GAUSS, x0-0.5*dx, x0+0.5*dx, g_gps, g_gws, l_gps_x, l_gws_x);
        // integrate H(distance) to determine VOF value
        for(int kk=0; kk<ORDER_GAUSS; kk++){
          for(int jj=0; jj<ORDER_GAUSS; jj++){
            for(int ii=0; ii<ORDER_GAUSS; ii++){
              double dist = l_gps_y[jj]-cy;
              VOF(i, j, k) +=
                l_gws_x[ii]*l_gws_y[jj]*l_gws_z[kk]
                *0.5*(1.+tanh(VOFBETA*dist/delta))
                /dx/dy/dz;
            }
          }
        }
      }
    }
  }
  common_free(g_gps);
  common_free(g_gws);
  common_free(l_gps_x);
  common_free(l_gws_x);
  common_free(l_gps_y);
  common_free(l_gws_y);
  common_free(l_gps_z);
  common_free(l_gws_z);
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        VOF(i, j, k) = fmin(1., VOF(i, j, k));
        VOF(i, j, k) = fmax(0., VOF(i, j, k));
      }
    }
  }
  interface_update_boundaries_vof(domain, vof);
  return 0;
}

static int load(const domain_t *domain, interface_t *interface){
  const int glisize = domain->glsizes[0];
  const int gljsize = domain->glsizes[1];
  const int glksize = domain->glsizes[2];
  const int   isize = domain->mysizes[0];
  const int   jsize = domain->mysizes[1];
  const int   ksize = domain->mysizes[2];
  const int ioffset = domain->offsets[0];
  const int joffset = domain->offsets[1];
  const int koffset = domain->offsets[2];
  // vof: [0:isize+1] x [1:jsize] x [1:ksize]
  double * restrict vof = interface->vof;
  const int glsizes[NDIMS] = {glksize, gljsize, glisize+2};
  const int mysizes[NDIMS] = {  ksize,   jsize,   isize+2};
  const int offsets[NDIMS] = {koffset, joffset, ioffset  };
  double *buf = common_calloc(mysizes[0] * mysizes[1] * mysizes[2], sizeof(double));
  const char *dirname = config.get_string("restart_dir");
  fileio_r_nd_parallel(dirname, "vof", NDIMS, glsizes, mysizes, offsets, buf);
  for(int cnt = 0, k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 0; i <= isize+1; i++){
        VOF(i, j, k) = buf[cnt];
        cnt++;
      }
    }
  }
  common_free(buf);
  interface_update_boundaries_vof(domain, vof);
  return 0;
}

#endif // NDIMS

interface_t *interface_init(const domain_t *domain){
  /* memory allocation */
  interface_t *interface = allocate(domain);
  /* init gaussian quadrature (including interval conversion */
  double *gps = interface->gps;
  double *gws = interface->gws;
  // compute roots and weights of gaussian quadrature, whose interval is [-1, 1]
  double max_residual = init_gaussian_quadrature(ORDER_GAUSS, gps, gws);
  // integral range is [-0.5, 0.5] instead of [-1, 1]
  convert_gaussian_quadrature(ORDER_GAUSS, -0.5, 0.5, gps, gws, gps, gws);
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if(myrank == 0){
    printf("Maximum residual of Gaussian quadrature: %.1e\n", max_residual);
  }
  /* initialise vof field */
  const bool restart_sim = config.get_bool("restart_sim");
  if(restart_sim){
    load(domain, interface);
  }else{
    init(domain, interface);
  }
  return interface;
}

