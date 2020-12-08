/* See intsort.h */
/* Last edited on 2005-06-05 13:33:38 by stolfi */

#include <intsort.h>

#include <stdlib.h>
#include <affirm.h>

void isrt_binssort(int *h, int n, int cmp(int x, int y), int sgn)
  { 
    int i;
    for (i = 1; i < n; i++)
      { /* Binary search to locate position {h[j]} of {h[i]} in {h[0..i-1]}: */
        int v = h[i];
        int r = 0, s = i-1;
        while (r <= s)
          { /* Now {h[0..r-1] <= h[i]}, {h[s+1..i-1] > h[i]}: */
            int m = (r + s)/2; 
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
