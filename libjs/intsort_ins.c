/* See intsort.h */
/* Last edited on 2023-03-18 11:23:16 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdlib.h>

#include <affirm.h>

#include <intsort.h>

void isrt_inssort(int32_t *h, int32_t n, int32_t cmp(int32_t x, int32_t y), int32_t sgn)
  { int32_t k;
    for (k = 1; k < n; k++)
      { int32_t v = h[k]; int32_t r = k;
        while ((r > 0) && (sgn*cmp(h[r-1],v) > 0)) { h[r] = h[r-1]; r--; }
        h[r] = v;
      }
  }
