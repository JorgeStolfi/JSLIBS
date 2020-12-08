/* See {hrn_more.h} */
/* Last edited on 2007-10-14 15:32:57 by stolfi */

#include <hrn_more.h>

void hrn_regular_simplex(int n, double p[], double h[])
  { double N = (double)n;
    double N1 = (double)(n+1);
    double SN1 = sqrt(N1);
    double b = 1/SN1;
    double m = (1 + b)/N;
    int i, j;
    int n1 = n+1;
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
          { int n1i = i*n1;
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
