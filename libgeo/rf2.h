/* rf2.h --- operations on points and vectors of R^2 (single-precision version) */
/* Last edited on 2024-11-15 00:34:20 by stolfi */

#ifndef rf2_H
#define rf2_H

#define _GNU_SOURCE
#include <stdint.h>

typedef struct rf2_t { float c[2]; } rf2_t;

rf2_t rf2_scale (double s, rf2_t* const a);
  /* Returns the vector {s * a}. */
  
rf2_t rf2_throw_cube(void);
  /* Returns a random vector in the cube {[-1 _ +1]^3} */

#endif
