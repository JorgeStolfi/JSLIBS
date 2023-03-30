/* See intsort.h */
/* Last edited on 2023-03-18 11:25:33 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdlib.h>

#include <affirm.h>

#include <intsort.h>

void isrt_binssort(int32_t *h, int32_t n, int32_t cmp(int32_t x, int32_t y), int32_t sgn)
  { 
    int32_t i;
    for (i = 1; i < n; i++)
      { /* Binary search to locate position {h[j]} of {h[i]} in {h[0..i-1]}: */
        int32_t v = h[i];
        int32_t r = 0, s = i-1;
        while (r <= s)
          { /* Now {h[0..r-1] <= h[i]}, {h[s+1..i-1] > h[i]}: */
            int32_t m = (r + s)/2; 
            if (sgn*cmp(v,h[m]) < 0)
              { s = m-1; }
            else
              { r = m+1; }
          }
        /* Displace successors and insert: */
        for (s = i; s > r; s--)
          { h[s] = h[s-1]; } 
        h[r] = v;
      }
  }
