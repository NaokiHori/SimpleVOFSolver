#if !defined(INTERFACE_H)
#define INTERFACE_H

#include "array.h"

// data type for N-dimensional vectors
typedef double vector_t[NDIMS];
// data type for (N+1)-dimensional vectors
//   to store additional element (segment),
//   which is stored at the last
typedef double normal_t[NDIMS + 1];

typedef struct {
  array_t vof;
  array_t ifrcx;
  array_t ifrcy;
  array_t dvof;
  array_t normal;
  array_t curv;
  array_t flxx;
  array_t flxy;
  array_t src[2];
  double tension;
} interface_t;

#endif // INTERFACE_H
