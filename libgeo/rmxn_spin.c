/* See rmxn_spin.h. */
/* Last edited on 2024-11-23 18:55:07 by stolfi */

#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include <rn.h>
#include <rmxn.h>
#include <rmxn_throw.h>
#include <affirm.h>

#include <rmxn_spin.h>
void rmxn_spin_rows(uint32_t m, uint32_t n, double A[], double M[])
  { /* Generate a random orthonormal {n×n} matrix {N}: */
    double N[n*n];
    rmxn_throw_ortho(n, N);
    /* Map each row of {A} by {N} (beware of aliasing between {A} and {M}): */
    double v[n];
    for (uint32_t i = 0;  i < m; i++)
      { rmxn_map_row(n, n, &(A[i*n]), N, v);
        rn_copy(n, v, &(M[i*n]));
      }
  }

void rmxn_spin_cols(uint32_t m, uint32_t n, double A[], double M[])
  { /* Generate a random orthonormal {m×m} matrix {N}: */
    double N[m*m];
    rmxn_throw_ortho(m, N);
    /* Map each col of {A} by {N} (beware of aliasing between {A} and {M}): */
    double a[m], v[m];
    for (uint32_t j = 0; j < n; j++)
      { rmxn_get_col(m, n, A, j, a);
        rmxn_map_col(m, m, a, N, v);
        rmxn_set_col(m, n, M, j, v);
      }
  }
