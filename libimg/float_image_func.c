/* See {float_image_func.h}. */
/* Last edited on 2018-03-04 22:42:49 by stolfilocal */

#include <math.h>
#include <limits.h>
#include <stdint.h>
#include <assert.h>
 
#include <bool.h>
#include <jsmath.h>
#include <affirm.h>
#include <float_image.h>
#include <float_image_func.h>

void float_image_func_get_index_range(double z, double r, int32_t n, int32_t *iLo, int32_t *iHi)
  {
    /* Compute float range with extra 0.5 margin (exclusive): */
    double zLo = fmax(-1000000.0, z - r - 0.4999999);
    double zHi = fmin(+1000000.0, z + r + 0.4999999);
    /* Convert range to integer: */
    int32_t lo = (int32_t)floor(zLo); 
    int32_t hi = (int32_t)floor(zHi); 
    if (n > 0) 
      { /* Clip to {[0..n-1]}: */
        if (lo < 0) { lo = 0; }
        if (hi >= n) { hi = n - 1; }
      }
    (*iLo) = lo; (*iHi) = hi;
  }

