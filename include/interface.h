#if !defined(INTERFACE_H)
#define INTERFACE_H

#include "domain.h"
#include "structure.h"

#define VOFMIN 1.e-8
#define VOFBETA 2.
#define ORDER_GAUSS 2

#if NDIMS == 2

typedef struct {
  double x;
  double y;
} dvof_t;

typedef struct {
  double a00;
  double a10;
  double a01;
} nrml_t;

struct interface_t_ {
  double *vof;
  dvof_t *dvof;
  nrml_t *normal;
  double *curv;
  double *voffluxx, *voffluxy;
  double *voffrcx, *voffrcy;
  double *vofsrca, *vofsrcb;
  double *gps, *gws;
};

#else // NDIMS == 3

typedef struct {
  double x;
  double y;
  double z;
} dvof_t;

typedef struct {
  double a000;
  double a100;
  double a010;
  double a001;
} nrml_t;

struct interface_t_ {
  double *vof;
  dvof_t *dvof;
  nrml_t *normal;
  double *curv;
  double *voffluxx, *voffluxy, *voffluxz;
  double *voffrcx, *voffrcy, *voffrcz;
  double *vofsrca, *vofsrcb;
  double *gps, *gws;
};

#endif // NDIMS

/* initialiser and finaliser */
extern interface_t *interface_init(const domain_t *domain);
extern int interface_finalise(interface_t *interface);

/* curvature tensor (normal and curvature) */
extern int interface_compute_curvature_tensor(const domain_t *domain, interface_t *interface);
/* compute surface force, csf model */
extern int interface_compute_surface_force(const domain_t *domain, interface_t *interface);
/* udpate */
extern int interface_update_vof(const domain_t *domain, const int rkstep, const double dt, const fluid_t *fluid, interface_t *interface);

#endif // INTERFACE_H
