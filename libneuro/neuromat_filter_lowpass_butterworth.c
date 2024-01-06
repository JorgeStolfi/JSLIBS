/* See {neuromat_filter_lowpass_butterworth.h}. */
/* Last edited on 2024-01-05 17:37:20 by stolfi */

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
    demand(isfinite(fc) && (fs > 0) && (fc != fs), "invalid {fc}");
    demand((gc > 0) && (gc < 1), "invalid {gc}");
    double rc = 1.0/(gc*gc) - 1.0;
    demand(rc > 0, "order would be infinite");
    if (rc <= 1) { return 1; }
    double tc = fabs(log(fc) - log(fs));
    assert(tc > 0);
    double oc = log(rc)/tc;
    assert(oc > 0);
    demand(oc <= (1 << 20), "order would be too big");
    int32_t ord = (int32_t)ceil(oc);
    return ord;
  }

double neuromat_filter_lowpass_butterworth_compute_fsup(double fs, int32_t ord)
  { demand(isfinite(fs) && (fs >= 0), "invalid {fs}");
    demand(ord >= 1, "invalid {ord}");
    if (fs == 0) 
      { return 0.0; }
    else
      { return fs*exp(36.0/ord); }
  }

double neuromat_filter_lowpass_butterworth_eval(double f, double fs, int32_t ord)
  {
    demand(fs >= 0, "invalid {fs}");
    demand(ord >= 1, "invalid {ord}");
    f = fabs(f);
    if (fs == 0) { return 0; } 
    if (f == 0) { return 1; }
    double z = 2.0*ord*(log(f) - log(fs));
    double w = NAN;
    if (z >= 71.0)
      { w = 0; }
    else 
      { if (z <= -35)
          { w =  1; }
        else if (z >= 37)
          { w = exp(-0.5*z); }
        else
          { w = sqrt(1/(1 + exp(z))); }
      }
    return w;
  }
  
