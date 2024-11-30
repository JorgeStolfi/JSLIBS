/* See gauss_elim_residual.h */
/* Last edited on 2024-11-25 01:50:54 by stolfi */

#include <stdint.h>
#include <assert.h>

#include <gauss_elim_residual.h>

void gauss_elim_residual
  ( uint32_t m,
    uint32_t n,
    double A[],
    uint32_t p,
    double B[],
    double X[],
    double R[]
  )
  { /* Compute {A X - B} with care: */
    
    uint32_t in0 = 0; /* {== i*n}. */
    uint32_t ipj = 0; /* {== i*p+j}. */
    for (uint32_t i = 0;  i < m; i++)
      { for (uint32_t j = 0;  j < p; j++)
          { double sum = 0.0, corr = 0.0;
            uint32_t kpj = j; /* {== k*p+j}. */
            for (uint32_t k = 0;  k <= n; k++)
              { double term = (k == n ? -B[ipj] : A[in0+k]*X[kpj]);
                /* Kahan's summation: */
                double tcorr = term - corr;
                double newSum = sum + tcorr;
                corr = (newSum - sum) - tcorr;
                sum = newSum;
                kpj += p;
              }
            R[ipj] = sum; ipj++;
          }
        in0 += n;
      }
  }
