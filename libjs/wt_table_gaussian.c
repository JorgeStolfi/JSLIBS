/* See wt_table_gaussian.h */
/* Last edited on 2023-11-25 12:11:18 by stolfi */

#define wt_table_gaussian_C_COPYRIGHT \
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

#include <wt_table_gaussian.h>

/* IMPLEMENTATIONS */

double_vec_t wt_table_gaussian_make(int32_t n, double sigma, double maxLoss)
  { if (n == 0)
      { /* Find {r} so that the omitted weight is at most {maxLoss}: */
        int32_t r = 0;
        while (wt_table_gaussian_loss(2*r + 1, sigma) > maxLoss) { r++; }
        assert(r > 0);
        n = 2*r + 1; 
      }
    /* Allocate and fill the table: */
    double_vec_t wt = double_vec_new(n);
    wt_table_gaussian_fill(wt.ne, sigma, wt.e, NULL);
    return wt;
  }
  
double wt_table_gaussian_loss(int32_t n, double sigma)
  { double r1 = 0.5*n, t1 = r1/sigma; 
    double ws = erfc(t1/M_SQRT2); 
    return ws;
  }
  
void wt_table_gaussian_fill(int32_t n, double sigma, double wt[], int32_t *stride_P)
  { demand(n > 0, "invalid table length");
    demand(sigma >= 0, "invalid {sigma}");
    int32_t stride;
    int32_t r = n/2;
    if ((sigma == 0) || (n <= 2))
      { for (int32_t k = 0; k < n; k++) { wt[k] = 0; }
        if ((n & 1) == 1)
          { wt[r] = 1.0; 
            stride = 1;
          }
        else
          { assert(n >= 2);
            wt[r-1] = 0.5;
            wt[r] = 0.5;
            stride = 2;
          }
       }
    else
      { assert ((sigma > 0) && (n > 2));
        /* Compute entries {wt[k]} and their sum {sumw}: */
        for (int32_t k = 0, j = n-1; k < n; k++, j--) 
          { wt[k] = (k <= j ? wt_table_gaussian_entry(n, k, sigma) : wt[j]); }
        stride = ((n == 3) && (wt[1] == 2*wt[0]) ? 1 : 0);
      }
    if (stride_P != NULL) { (*stride_P) = stride; }
  }
  
double wt_table_gaussian_entry(int32_t n, int32_t k, double sigma)
  { /* Integral of Gaussian within interval {k} of {n} intervals centered at origin: */
    double r0 = ((double)k) - 0.5*n;
    double r1 = r0 + 1.0; 
    double wr = gauss_distr_integral(r0, r1, 0.0, sigma);
    return wr;
  }
