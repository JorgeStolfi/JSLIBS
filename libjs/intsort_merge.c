/* See intsort.h */
/* Last edited on 2024-11-17 15:50:53 by stolfi */

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <affirm.h>
#include <intmerge.h>

#include <intsort.h>

void isrt_mergesort (int32_t *h, uint32_t n, int32_t cmp(int32_t x, int32_t y), int32_t sgn)
  {
    auto void sort(int32_t *a, int32_t *b);
    
    if (n >= 2) { sort(h, h+n); }
    return;
    
    void sort(int32_t *a, int32_t *b)
      { assert(a < b);
        uint32_t nab = (uint32_t)(b - a);
        if (nab <= 1) 
          { /* Nothing to do. */ }
        else
          { int32_t *m = a + nab/2;
            if (a < m) { sort (a, m); }
            if (m < b) { sort (m, b); }
            imrg_merge(a, m, b, cmp, sgn);
          }
      }
  }

/* Created by J. Stolfi, Unicamp, Oct/2004. */
