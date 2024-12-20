/* See sym_eigen_tridiag.h */
/* Last edited on 2024-12-05 21:42:02 by stolfi */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>

#include <sym_eigen_tred1.h>
#include <sym_eigen_tred2.h>

#include <sym_eigen_tridiag.h>

#define Pr fprintf
#define Er stderr
  /* Shoter names. */

void sym_eigen_tridiag
  ( uint32_t n,
    double A[],
    double d[],
    double e[],
    double R[]
  )
  {
    /* The algorithm is based on the EISPACK routines "tred1.f" and
      "tred2.f" by Burton S. Garbow, Argonne National Laboratory
      (aug/1983); originally from the ALGOL procedures {tred1} and
      {tred2}, Num. Math. 11, 181-195(1968) by Martin, Reinsch, and
      Wilkinson. See Handbook for Auto. Comp., vol.II-Linear Algebra,
      212-226(1971). Re-implemented in C by Jorge Stolfi, Unicamp
      (dec/2002). */
      
    if (n == 0) 
      { /* Nothing to do: */
      }
    else if (n == 1)
      { /* Trivial case: */
        d[0] = A[0];
        e[0] = 0;
        if (R != NULL) { R[0] = 1.0; }
      }
    else if (n == 2)
      { /* {A} is already tridiagonal, so {R} is identity: */
        d[0] = A[0]; d[1] = A[3];
        e[0] = 0; e[1] = A[2];
       if (R != NULL) {  R[0] = R[3] = 1.0; R[1] = R[2] = 0.0; }
      }
    else
      { if (R == NULL)
          { sym_eigen_tred1(n, A, d, e, e); }
        else
          { sym_eigen_tred2(n, A, d, e, R); }
      }
  }
