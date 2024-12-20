/* See sym_eigen.h */
/* Last edited on 2024-12-05 21:45:44 by stolfi */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>

#include <sym_eigen_tridiag.h>
#include <sym_eigen_tql1.h>
#include <sym_eigen_tql2.h>

#include <sym_eigen.h>

void sym_eigen
  ( uint32_t n,
    double A[],
    double d[],
    double R[],
    uint32_t *nev_P
  )
  {
    if (n == 0) 
      { /* Nothing to do: */ 
        (*nev_P) = 0;
        return;
      }
    else if (n == 1)
      { /* Trivial case: */
        d[0] = A[0];
        if (R != NULL) { R[0] = 1.0; }
        (*nev_P) = 1;
        return;
      }
    else
      { double e[n];
        sym_eigen_tridiag(n, A, d, e, R);
        if (R == NULL)
          { sym_eigen_tql1(n, d, e, nev_P); }
        else
          { sym_eigen_tql2(n, d, e, R, nev_P); }
      }
  }
