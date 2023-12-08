/* See wt_table_triangular.h */
/* Last edited on 2023-11-26 05:47:57 by stolfi */

#define wt_table_triangular_C_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <values.h>

#include <affirm.h>
#include <gauss_distr.h>

#include <wt_table_triangular.h>

double_vec_t wt_table_triangular_make(int32_t n)
  { /* Allocate and fill the table: */
    double_vec_t wt = double_vec_new(n);
    wt_table_triangular_fill(n, wt.e, NULL);
    return wt;
  }
  
void wt_table_triangular_fill(int32_t n, double wt[], int32_t *stride_P)
  { /* Build triangular table and compute its sum: */
    int32_t m = (n-1)/2;
    for (int32_t k = 0; k < n; k++)
      { wt[k] = (k <= m ? k + 1.0 : wt[n-1-k]); }
    if (stride_P != NULL) { (*stride_P) = (n+1)/2; }
  }
