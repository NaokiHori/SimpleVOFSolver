#if !defined(SAVE_H)
#define SAVE_H

#include "domain.h"
#include "fluid.h"
#include "temperature.h"
#include "interface.h"

// next time to trigger output
extern double save_next;

extern int save(const domain_t *domain, const int step, const double time, const fluid_t *fluid, const temperature_t *temperature, const interface_t *interface);

#endif // SAVE_H
