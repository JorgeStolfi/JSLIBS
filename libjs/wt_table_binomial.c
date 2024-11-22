/* See wt_table_binomial.h */
/* Last edited on 2024-11-16 10:25:46 by stolfi */

#define wt_table_binomial_C_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <values.h>

#include <affirm.h>
#include <gauss_distr.h>

#include <wt_table_binomial.h>

double_vec_t wt_table_binomial_make(uint32_t n)
  { /* Allocate and fill the table: */
    double_vec_t wt = double_vec_new(n);
    return wt;
  }

void wt_table_binomial_fill(uint32_t n, double wt[], uint32_t *stride_P)
  { /* Repeated convolutions (Pascal's triangle): */
    demand(n >= 1, "invalid table size {n}");
    wt[0] = 1;
    for (int32_t k = 1; k < n; k++)
      { wt[k] = wt[k-1];
        for (int32_t j = k-1; j > 0; j--) { wt[j] = wt[j] + wt[j-1]; }
      }
    if (stride_P != NULL) { (*stride_P) = (n == 1 ? 1 : 2); }
  }
