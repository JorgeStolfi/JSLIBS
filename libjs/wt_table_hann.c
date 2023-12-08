/* See wt_table_hann.h */
/* Last edited on 2023-11-25 12:13:19 by stolfi */

#define wt_table_hann_C_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <values.h>

#include <affirm.h>
#include <gauss_distr.h>

#include <wt_table_hann.h>
   
void wt_table_hann_fill(int32_t n, double flat, double wt[], int32_t *stride_P)
  { demand((flat >= 0) && (flat <= 1.0), "invalid {flat} fraction");
    /* Compute parameters: */
    double c = 0.5*(1-flat)*(n-1);  /* index (fractional) where the weight becomes 1. */
    for (int32_t k = 0; k < n; k++)
      { if (k > n-1-k)
          { wt[k] = wt[n-1-k]; }
        else if (k < c)
          { wt[k] = 0.5*(1 + cos(M_PI*(k-c)/(c+1))); }
        else
          { wt[k] = 1.0; }
      }
    if (stride_P != NULL)
      { int32_t r = (int32_t)floor(c + 1.0e-13);
        (*stride_P) = (fabs(c - (double)r) < 1.0e-12 ? n - r : 0);
      }
  }
   
double_vec_t wt_table_hann_make(int32_t n, double flat)
  { /* Allocate and fill the table: */
    double_vec_t wt = double_vec_new(n);
    wt_table_hann_fill(n, flat, wt.e, NULL);
    return wt;
  }
