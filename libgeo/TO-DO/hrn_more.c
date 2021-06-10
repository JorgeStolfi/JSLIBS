/* See {hrn_more.h} */
/* Last edited on 2021-06-09 20:13:49 by jstolfi */

#include <stdint.h>
#include <hrn_more.h>

void hrn_regular_simplex(int32_t n, double p[], double h[])
  { double N = (double)n;
    double N1 = (double)(n+1);
    double SN1 = sqrt(N1);
    double b = 1/SN1;
    double m = (1 + b)/N;
    int32_t i, j;
    int32_t n1 = n+1;
    if (p != NULL)
      { /* Set the matrix {p}: */
        /* ... */
      }
    if (h != NULL)
      { /* Set the matrix {h}: */
        double r = b;
        double s = 1;
        double u = ...;
        double v = ...;
        double w = ...;
        for (i = 0; i <= n; i++) 
          { int32_t n1i = i*n1;
            if (i == 0)
              { /* Set the first row to {<r,s,s,...,s>}: */
                h[n1i] = r; /* Weight. */
                for (j = 1; j <= n; j++) { h[n1i + j] = s; }
              }
            else
              { /* Set row {i} to {<u,v,v,...,w,...v>}: */
                h[n1i] = u; /* Weight. */
                for (j = 1; j <= n; j++) { h[n1i + j] = (i == j ? w : v); }
              }
          }
      }
   }
