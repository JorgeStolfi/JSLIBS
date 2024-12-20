/* See {rmxn_transform_quadratic.h}. */
/* Last edited on 2024-12-05 15:54:13 by stolfi */

#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
 
#include <bool.h>
#include <rmxn.h>
#include <jsmath.h>
#include <sym_eigen.h>
#include <sym_eigen_test_tools.h>
#include <affirm.h>

#include <rmxn_transform_quadratic.h>

void rmxn_transform_quadratic
  ( uint32_t n,
    double E[],
    double e[],
    uint32_t m,
    double U[],
    double F[],
    double f[]
  )
  { 
    bool_t debug = TRUE;
    if (debug) { fprintf(stderr, "  > %s\n", __FUNCTION__); }
    demand(n >= 0, "invalid {n}");
    demand(m >= 0, "invalid {m}");
    if (m == 0) { return; }

    if (debug) { fprintf(stderr, "    ... Computing the matrix {H = U E} ...\n"); }
    double H[m*n];
    rmxn_mul(m, n, n, U, E, H);

    if (debug) { fprintf(stderr, "    ... computing the metric matrix {M} for {\\QF} ...\n"); }
    double M[m*m];
    for (uint32_t r = 0;  r < m; r++)
      { double *hr = &(H[r*n]);
        for (uint32_t s = 0;  s <= r; s++)
          { double *hs = &(H[s*n]);
            double sum = 0.0;
            for (uint32_t j = 0; j < n; j++)
              { double hrj = hr[j];
                double hsj = hs[j];
                double ej = e[j];
                sum += hrj*ej*hsj; 
              }
            M[r*m + s] = sum;
            M[s*m + r] = sum; /* Diag is assigned twice, but OK. */
          }
      }
    
    if (debug) { fprintf(stderr, "    ... computing the eigen decomp {F,f} of {M} ...\n"); }
    uint32_t nev; /* Number of eigenvalues successfully found. */
    sym_eigen(m, M, f, F, &nev);
    if (debug)
      { double tiny = 1.0e-14*rmxn_norm(m,m,M)/m;
        fprintf(stderr, "    nev = %d (should be %d)\n", nev, m);
        sym_eigen_test_tools_print_matrix(stderr, "%18.12f", m, m,m, M, tiny);
      }
    demand(nev == m, "failed to determine eigenvalues of {M}");
    if (debug) { fprintf(stderr, "  < %s\n", __FUNCTION__); }
  }

