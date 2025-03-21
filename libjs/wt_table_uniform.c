/* See wt_table_uniform.h */
/* Last edited on 2024-11-16 10:31:44 by stolfi */

#define wt_table_uniform_C_COPYRIGHT \
  "Copyright � 2006  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <values.h>

#include <affirm.h>
#include <gauss_distr.h>

#include <wt_table_uniform.h>

double_vec_t wt_table_uniform_make(uint32_t n, double val)
  { /* Allocate and fill the table: */
    double_vec_t wt = double_vec_new(n);
    wt_table_uniform_fill(wt.ne, val, wt.e, NULL);
    return wt;
  }

void wt_table_uniform_fill(uint32_t n, double val, double wt[], uint32_t *stride_P)
  { demand(n >= 1, "invalid table size {n}");
    demand(val > 0, "invalid table value {val}");
    for (uint32_t k = 0;  k < n; k++) { wt[k] = val; }
    if (stride_P != NULL) { (*stride_P) = n; }
  }
