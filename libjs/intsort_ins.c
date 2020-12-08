/* See intsort.h */
/* Last edited on 2005-06-05 13:33:02 by stolfi */

#include <intsort.h>

#include <stdlib.h>
#include <affirm.h>

void isrt_inssort(int *h, int n, int cmp(int x, int y), int sgn)
  { int k;
    for (k = 1; k < n; k++)
      { int v = h[k]; int r = k;
        while ((r > 0) && (sgn*cmp(h[r-1],v) > 0)) { h[r] = h[r-1]; r--; }
        h[r] = v;
      }
  }
