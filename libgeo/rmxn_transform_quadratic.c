/* See {rmxn_transform_quadratic.h}. */
/* Last edited on 2022-01-03 17:10:20 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
 
#include <bool.h>
#include <rmxn.h>
#include <jsmath.h>
#include <sym_eigen.h>
#include <affirm.h>

#include <rmxn_transform_quadratic.h>

void rmxn_transform_quadratic(int32_t n, double E[], double e[], int32_t m, double U[], double F[], double f[])
  { 
    bool_t debug = TRUE;
    
    demand(n >= 0, "invalid {n}");
    demand(m >= 0, "invalid {m}");
    if (m == 0) { return; }

    if (debug) { fprintf(stderr, "... Computing the matrix {H = U E} ...\n"); }
    double H[m*n];
    rmxn_mul(m, n, n, U, E, H);

    if (debug) { fprintf(stderr, "... computing the metric matrix {M} for {\\EF} ...\n"); }
    double M[m*m];
    for (int32_t r = 0; r < m; r++)
      { double *hr = &(H[r*n]);
        for (int32_t s = 0; s <= r; s++)
          { double *hs = &(H[s*n]);
            double sum = 0.0;
            for (int32_t i = 0; i < n; i++)
              { double hrij = hr[i];
                double hsij = hs[i];
                double ei = e[i];
                sum += hrij*hsij*ei; 
              }
            M[r*m + s] = sum;
            M[s*m + r] = sum; /* Diag is assigned twice, but OK. */
          }
      }
    
    if (debug) { fprintf(stderr, "... computing the eigen decomp {F,f} of {M} ...\n"); }
    /* Convert {M} to tridiag with diagonal {d[0..m-1]} and subdiagonal {t[1..m-1]}: */
    double d[m]; /* Diagonal elements of tridiagonal matrix. */
    double t[m]; /* Sub-diagonal elements of temporary tridiagonal matrix. */
    syei_tridiagonalize(m, M, d, t, F);
    /* Compute eigenvalues and eigenvectors from {F} and tridiag matrix: */
    int32_t p; /* Number of eigenvalues computed. */
    int32_t absrt = 0; /* Sort eigenvalues by signed value. */
    syei_trid_eigen(m, d, t, F, &p, absrt);
    /* Check that all eigenvalues were computed: */
    demand(p == m, "failed to determine eigenvalues of {M}");
    /* Copy the eigenvalues {d[0..m-1]} to {f[0..m-1]}: */
    for (int32_t k = 0; k < m; k++) { f[k] = d[k]; } 
    
  }

