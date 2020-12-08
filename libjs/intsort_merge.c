/* See intsort.h */
/* Last edited on 2013-10-25 01:36:13 by stolfilocal */

#include <intsort.h>

#include <stdlib.h>
#include <intmerge.h>
#include <affirm.h>

void isrt_mergesort (int *h, int n, int cmp(int x, int y), int sgn)
  {
    auto void sort(int *a, int *b);

    void sort(int *a, int *b)
    { int nab = (int)(b - a);
        if (nab <= 1) 
          { /* Nothing to do. */ }
        else
          { int *m = a + nab/2;
            if (a < m) { sort (a, m); }
            if (m < b) { sort (m, b); }
            imrg_merge(a, m, b, cmp, sgn);
          }
      }
    sort(h, h+n);
  }

/* Created by J. Stolfi, Unicamp, Oct/2004. */
