/* See intsort.h */
/* Last edited on 2024-11-23 06:07:31 by stolfi */

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <affirm.h>

/* FOR DEBUGGING: */
#include <stdio.h>
#include <bool.h>

#include <intsort.h>

void isrt_binssort(int32_t *h, uint32_t n, int32_t cmp(int32_t x, int32_t y), int32_t sgn)
  { bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "  > %s n = %u ---\n", __FUNCTION__, n); }
    for (int32_t i = 1;  i < n; i++)
      { /* Binary search to locate position {h[j]} of {h[i]} in {h[0..i-1]}: */
        if (debug) { fprintf(stderr, "   i = %u\n", i); }
        int32_t v = h[i];
        int32_t r = 0, s = i;
        while (r < s)
          { /* Now {h[0..r-1] <= h[i]}, {h[s..i-1] > h[i]}. */
            /* Locate {v} among {h[r..s-1]}: */
            int32_t m = (r + s)/2; 
            if (debug) { fprintf(stderr, "     r = %u m = %u s = %u\n", r, m, s); }
            assert((r <= m) && (m < s));
            if (sgn*cmp(v,h[m]) < 0)
              { s = m; }
            else
              { r = m+1; }
          }
        /* Now {s == r} and  {h[0..r-1] <= h[i]}, {h[r..i-1] > h[i]}. */
        /* Displace {h[r..i-1]} and insert: */
        for (s = i; s > r; s--) { h[s] = h[s-1]; } 
        h[r] = v;
      }
    if (debug) { fprintf(stderr, "  < %s ---\n", __FUNCTION__); }
  }
