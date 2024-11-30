/* See gauss_elim_determinant.h */
/* Last edited on 2024-11-25 01:43:30 by stolfi */

#include <stdint.h>
#include <assert.h>

#include <bool.h>

#include <gauss_elim_triangularize.h>
#include <gauss_elim_determinant.h>

double gauss_elim_determinant(uint32_t m, uint32_t n, double A[], uint32_t q)
  { if ((q > m) || (q > n)) { return 0.0; }
    /* Make a work copy of the first {q} rows and columns of {A}: */
    double M[q*q];
    
    for (int32_t i = 0; i < q; i ++)
      { for (int32_t j = 0;  j < q; j++) { M[i*(int32_t)q + j] = A[i*(int32_t)n + j]; } }
    
    /* Triangularize and get determinant: */
    gauss_elim_triangularize(q, q, M, FALSE, 0.0);
    return gauss_elim_determinant_triang(q, q, M, q);
  }

double gauss_elim_determinant_triang(uint32_t m, uint32_t n, double M[], uint32_t q)
  { if ((q > m) || (q > n)) { return 0.0; }
    double det = 1.0;
    for (int32_t i = 0; (i < q) && (det != 0); i++)
      { det *= M[i*(int32_t)n + i];}
    return det;
  }
