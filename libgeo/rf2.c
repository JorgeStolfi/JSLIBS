/* See rf2.h. */
/* Last edited on 2024-11-20 13:01:40 by stolfi */

#include <stdint.h>
#include <assert.h>

#include <jsrandom.h>

#include <rf2.h>

#define N 2

rf2_t rf2_scale (double s, rf2_t* const a)
  { rf2_t r;
    r.c[0] = (float)(s * a->c[0]);
    r.c[1] = (float)(s * a->c[1]);
    return r;
  }

rf2_t rf2_throw_cube (void)
  { rf2_t r;
    for (uint32_t i = 0;  i < N; i++)
      { r.c[i] = (float)(2.0*drandom() - 1.0); }
    return r;
  }
