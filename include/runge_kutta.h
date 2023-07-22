#if !defined(RUNGE_KUTTA_H)
#define RUNGE_KUTTA_H

#include <stdint.h>

// Runge-Kutta configurations
// indices
extern const uint_fast8_t rk_a; // 0
extern const uint_fast8_t rk_b; // 1
extern const uint_fast8_t rk_g; // 2
// coefficients of three-step RK scheme
typedef double rkcoef_t[3];
extern const rkcoef_t rkcoefs[3];

#endif // RUNGE_KUTTA_H
