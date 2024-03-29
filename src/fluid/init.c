#include <math.h>
#include "param.h"
#include "memory.h"
#include "config.h"
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "fileio.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/fluid/p.h"
#include "array_macros/fluid/t.h"
#include "array_macros/fluid/psi.h"
#include "array_macros/fluid/srcux.h"
#include "array_macros/fluid/srcuy.h"
#include "array_macros/fluid/srct.h"

/**
 * @brief allocate members
 * @param[in]  domain : information about domain decomposition and size
 * @param[out] fluid  : structure storing flow fields and auxiliary buffers
 * @return            : error code
 */
static int allocate(
    const domain_t * domain,
    fluid_t * fluid
){
  // velocity
  if(0 != array.prepare(domain, UX_NADDS, sizeof(double), &fluid->ux )) return 1;
  if(0 != array.prepare(domain, UY_NADDS, sizeof(double), &fluid->uy )) return 1;
  // pressure and scalar potential
  if(0 != array.prepare(domain, P_NADDS,   sizeof(double), &fluid->p  )) return 1;
  if(0 != array.prepare(domain, PSI_NADDS, sizeof(double), &fluid->psi)) return 1;
  // temperature
  if(0 != array.prepare(domain, T_NADDS, sizeof(double), &fluid->t)) return 1;
  // Runge-Kutta source terms
  for(size_t n = 0; n < 3; n++){
    if(0 != array.prepare(domain, SRCUX_NADDS, sizeof(double), &fluid->srcux[n])) return 1;
    if(0 != array.prepare(domain, SRCUY_NADDS, sizeof(double), &fluid->srcuy[n])) return 1;
    if(0 != array.prepare(domain, SRCT_NADDS,  sizeof(double), &fluid->srct [n])) return 1;
  }
  return 0;
}

static void report(
    const sdecomp_info_t * info,
    const fluid_t * fluid
){
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(info, &myrank);
  if(root == myrank){
    printf("FLUID\n");
    printf("\tRa: % .7e\n", fluid->Ra);
    printf("\tPr: % .7e\n", fluid->Pr);
    printf("\tMomentum    diffusivity: % .7e\n", fluid->m_dif);
    printf("\tTemperature diffusivity: % .7e\n", fluid->t_dif);
    printf("\tdiffusive treatment in x: %s\n", param_m_implicit_x ? "implicit" : "explicit");
    printf("\tdiffusive treatment in y: %s\n", param_m_implicit_y ? "implicit" : "explicit");
    fflush(stdout);
  }
}

/**
 * @brief constructor of the structure
 * @param[in]  dirname_ic : name of directory in which initial flow fields are stored
 * @param[in]  domain     : information about domain decomposition and size
 * @param[out]            : structure being allocated and initalised
 * @return                : (success) 0
 *                          (failure) non-zero value
 */
int fluid_init(
    const char dirname_ic[],
    const domain_t * domain,
    fluid_t * fluid
){
  // allocate arrays
  if(0 != allocate(domain, fluid)) return 1;
  // load flow fields
  if(0 != array.load(domain, dirname_ic, "ux", fileio.npy_double, &fluid->ux)) return 1;
  if(0 != array.load(domain, dirname_ic, "uy", fileio.npy_double, &fluid->uy)) return 1;
  if(0 != array.load(domain, dirname_ic,  "p", fileio.npy_double, &fluid-> p)) return 1;
  if(0 != array.load(domain, dirname_ic,  "t", fileio.npy_double, &fluid-> t)) return 1;
  // impose boundary conditions and communicate halo cells
  fluid_update_boundaries_ux(domain, &fluid->ux);
  fluid_update_boundaries_uy(domain, &fluid->uy);
  fluid_update_boundaries_p(domain, &fluid->p);
  fluid_update_boundaries_t(domain, &fluid->t);
  // compute diffusivities
  if(0 != config.get_double("Pr", &fluid->Pr)) return 1;
  if(0 != config.get_double("Ra", &fluid->Ra)) return 1;
  fluid->m_dif = 1. * sqrt(fluid->Pr) / sqrt(fluid->Ra);
  fluid->t_dif = 1. / sqrt(fluid->Pr) / sqrt(fluid->Ra);
  report(domain->info, fluid);
  return 0;
}

