/* See {neuromat_filter_lowpass_gauss.h}. */
/* Last edited on 2023-12-15 22:54:23 by stolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <jsmath.h>

#include <neuromat_filter_lowpass_gauss.h>

double neuromat_filter_lowpass_gauss_compute_sigma(double fc, double gc)
  { demand(isfinite(gc) && (gc > 0) && (gc < 1), "invalid {gc}");
    double sigma = fc/sqrt(-2*log(gc));
    return sigma;
  }

double neuromat_filter_lowpass_gauss(double f, double sigma)
  { 
    demand(isfinite(sigma) && (sigma > 0), "invalid {sigma}");
    double zf = f/sigma;
    if (fabs(zf) > 9.0)
      { return 0.0; }
    else
      { return exp(-0.5*zf*zf); }
  } 
