/* See wt_table_quantize.h */
/* Last edited on 2023-11-04 18:43:42 by stolfi */

#define wt_table_quantize_C_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <values.h>

#include <affirm.h>
#include <jsmath.h>

#include <wt_table_quantize.h>

int32_t wt_table_quantize(int32_t n, double wd[], int32_t wi_min, int32_t wi_sum, int32_t wi[])
  { 
    demand(wi_min >= 0, "invalid negative {wi_min}");
    demand((wi_sum >= 0) && (wi_sum <= wt_table_quantize_WI_SUM_MAX), "invalid {wi_sum}");
    /* Compute the sum {wd_sum} of absolute weights: */
    double wd_sum = 0, wd_corr = 0;
    double wd_min = +INF;
    for (int32_t k = 0; k < n; k++)
      { double wdak = fabs(wd[k]);
        demand(isfinite(wdak), "weights cannot be infinite or NAN");
        if (wdak < wd_min) { wd_min = wdak; }
        /* Kahan's summation: */
        double tmp_corr = wdak - wd_corr;
        double tmp_sum = wd_sum + tmp_corr;
        wd_corr = (tmp_sum - wd_sum) - tmp_corr;
        wd_sum = tmp_sum;
      }
    if (wd_sum == 0)
      { /* All weights are zero: */
        for (int32_t k = 0; k < n; k++) { wi[k] = 0; }
        return 0;
      }

    if ((wi_sum == 0) && (wi_min == 0)) { wi_sum = wt_table_quantize_WI_SUM_MAX; }
    /* Choose tentative the scale factor, assuming worst-case rounding: */
    double scale = NAN;
    if (wi_sum == 0)
      { assert(wi_min != 0);
        scale = ((double)wi_min - 0.49999)/wd_min;
        /* But don't let {wi_sum} get too big: */
        double scale_max = ((double)wt_table_quantize_WI_SUM_MAX - n*wi_min)/(wd_sum + 0.5*n);
        demand(scale <= scale_max, "{wi_min} too large");
      }
    else
      { demand(wi_sum > n*wi_min, "{wi_min} too large for requested {wi_sum}");
        scale = ((double)wi_sum - n*wi_min)/(wd_sum + 0.5*n);
      }
    
    /* Round entries, compute actual sum of abs values: */
    int32_t wi_sum_cmp = 0;
    for (int32_t k = 0; k < n; k++) 
      { double wdk = wd[k];
        if (wdk != 0) 
          { int32_t wiak = (int32_t)floor(scale * fabs(wdk)); 
            if (wiak < wi_min) { wiak = wi_min; }
            wi[k] = (wdk < 0 ? -wiak : +wiak); 
            assert(wi_sum_cmp <= wt_table_quantize_WI_SUM_MAX - wiak); 
            wi_sum_cmp += wiak;
          }
        else
          { wi[k] = 0; }
      }
    assert(wi_sum_cmp != 0);

    return wi_sum_cmp;
  }
