#include "common.h"
#include "interface.h"


int interface_finalise(interface_t *interface){
  common_free(interface->vof);
  common_free(interface->dvof);
  common_free(interface->normal);
  common_free(interface->voffluxx);
  common_free(interface->voffluxy);
  common_free(interface->voffrcx);
  common_free(interface->voffrcy);
  common_free(interface->vofsrca);
  common_free(interface->vofsrcb);
  common_free(interface->gps);
  common_free(interface->gws);
  common_free(interface);
  return 0;
}

