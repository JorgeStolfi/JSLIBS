/* See {neuromat_filter_bandpass_log_gauss.h}. */
/* Last edited on 2023-12-15 10:42:05 by stolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <jsmath.h>

#include <neuromat_filter_bandpass_log_gauss.h>

double neuromat_filter_bandpass_log_gauss_compute_sigma(double fm, double fc, double gc)
  { demand(fm > 0, "invalid {fm}");
    demand((fc > 0) && (fc != fm), "invalid {fc}");
    demand((gc > 0) && (gc < 1), "invalid {gc}");
    double sigma = (log(fc) - log(fm))/sqrt(-2*log(gc));
    return sigma;
  }
  
double neuromat_filter_bandpass_log_gauss(double f, double fm, int32_t np, double sigma)
  { demand(fm > 0, "invalid {fm}");
    demand((sigma != 0), "invalid {sigma}");
    demand(f >= 0, "invalid {f}");
    demand(np >= 1, "invalid {np}");
    if (f == 0) { return 0; }
    double gtot = 0;
    for (int32_t kp = 0; kp < np; kp++)
      { double zf = fabs((log(f) - log(fm))/sigma - 2.0*kp);
        if (zf < 9.0) { double g = exp(-0.5*zf*zf); gtot += g; }
      }
    return gtot;
  } 

