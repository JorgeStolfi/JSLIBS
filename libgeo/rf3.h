/* rf3.h --- operations on points and vectors of R^3 (single-precision version) */
/* Last edited on 2024-12-05 10:28:19 by stolfi */

#ifndef rf3_H
#define rf3_H

#include <stdint.h>

typedef struct rf3_t { float c[3]; } rf3_t;

rf3_t rf3_add (rf3_t* const a, rf3_t* const b);
  /* Returns the vector {a + b}. */

rf3_t rf3_sub (rf3_t* const a, rf3_t* const b);
  /* Returns the vector {a - b}. */

rf3_t rf3_scale (double s, rf3_t* const a);
  /* Returns the vector {s * a}. */

rf3_t rf3_mix (double s, rf3_t* const a, double t, rf3_t* const b);
  /* Returns the vector {s * a + t * b}. */

float rf3_max (rf3_t* const a);
  /* Returns the maximum coordinate of {a}. */

rf3_t rf3_rot_axis (rf3_t* const a, uint32_t i, uint32_t j, double ang);
  /* Returns {a} after a rotation that moves Cartesian axis {i} towards 
    axis {j} by {ang} radians, leaving all other coordinates unchanged. */

rf3_t rf3_rot_gen (rf3_t* const a, rf3_t* const d, double ang);
  /* Returns {a} after a rotation by {ang} radians around the vector {r}
    in the sense of the right-hand rule. */

double rf3_norm (rf3_t* const a);
  /* Returns the Euclidean length of {a}. */
 
rf3_t rf3_throw_cube (void);
  /* Sets {r} to a uniformly random point of the 3-cube (square) {[-1 _ +1]^3}. */

#endif
