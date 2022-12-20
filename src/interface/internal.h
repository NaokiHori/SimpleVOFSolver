#if !defined(INTERFACE_INTERNAL_H)
#define INTERFACE_INTERNAL_H

#include "domain.h"
#include "structure.h"

/* Gaussian quadrature stuffs */
extern double init_gaussian_quadrature(const int n, const double *gps, const double *gws);
extern int convert_gaussian_quadrature(const int n, const double xm, const double xp, const double *gps_bef, const double *gws_bef, double *gps_aft, double *gws_aft);


/* common functions, surface function */
extern double psi(const double a10, const double a01, const double x, const double y);
extern double H(const double a00, const double a10, const double a01, const double x, const double y);


extern int interface_update_boundaries_vof(const domain_t *domain, double *vof);

#endif // INTERFACE_INTERNAL_H
