/* See intsort.h */
/* Last edited on 2023-03-18 11:22:44 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdlib.h>

#include <affirm.h>
#include <intmerge.h>

#include <intsort.h>

void isrt_mergesort (int32_t *h, int32_t n, int32_t cmp(int32_t x, int32_t y), int32_t sgn)
  {
    auto void sort(int32_t *a, int32_t *b);

    void sort(int32_t *a, int32_t *b)
    { int32_t nab = (int32_t)(b - a);
        if (nab <= 1) 
          { /* Nothing to do. */ }
        else
          { int32_t *m = a + nab/2;
            if (a < m) { sort (a, m); }
            if (m < b) { sort (m, b); }
            imrg_merge(a, m, b, cmp, sgn);
          }
      }
    sort(h, h+n);
  }

/* Created by J. Stolfi, Unicamp, Oct/2004. */
