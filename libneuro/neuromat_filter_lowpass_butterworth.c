/* See {neuromat_filter_lowpass_butterworth.h}. */
/* Last edited on 2023-12-16 00:59:24 by stolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <complex.h>

#include <affirm.h>
#include <bool.h>
#include <jsmath.h>

#include <neuromat_filter_lowpass_butterworth.h>

int32_t neuromat_filter_lowpass_butterworth_compute_order(double fs, double fc, double gc)
  { demand(isfinite(fs) && (fs > 0), "invalid {fs}");
    demand(isfinite(fc) && (fc > fs), "invalid {fc}");
    demand((gc > 0) && (gc < 1), "invalid {gc}");
    double rc = 1.0/(gc*gc) - 1.0;
    demand(rc > 0, "order would be infinite");
    if (rc <= 1) { return 1; }
    double tc = log(fc) - log(fs);
    assert(tc > 0);
    double oc = log(rc)/tc;
    assert(oc > 0);
    demand(oc <= (1 << 20), "order would be too big");
    int32_t ord = (int32_t)ceil(oc);
    return ord;
  }

double neuromat_filter_lowpass_butterworth(double f, double fs, int32_t ord)
  {
    demand(fs >= 0, "invalid {fs}");
    demand(f >= 0, "invalid {f}");
    if (f == 0) { return 1; }
    if (fs == 0) { return 0; } 
    double t = log(f) - log(fs);
    if (ord*t >= 74.0)
      { return 0; }
    else 
      { if (ord*t <= -37)
          { return 1; }
        else if (t >= 37)
          { return exp((0.5*ord)*t); }
        else
          { return sqrt(1/(1 + exp(ord*t))); }
      }
  }
  
