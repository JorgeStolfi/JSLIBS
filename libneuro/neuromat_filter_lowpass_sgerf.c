/* See {neuromat_filter_lowpass_sgerf.h}. */
/* Last edited on 2024-01-05 17:33:02 by stolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <jsmath.h>

#include <neuromat_filter_lowpass_sgerf.h>

double neuromat_filter_lowpass_sgerf_compute_sigma(double fm, double fc, double gc)
  { demand(isfinite(fm) && (fm > 0), "invalid {fm}");
    demand(isfinite(fc) && (fc > 0) && (fc != fm), "invalid {fc}");
    demand((gc > 0) && (gc < 1.0) && (gc != 0.5), "invalid {gc}");
    double den = fabs(erf_inv(1 - 2*gc));
    assert(! isnan(den));
    double sigma;
    if (den == +INF)
      { sigma = 0; }
    else if (den == 0)
      { sigma = +INF; }
    else
      { sigma = fabs(log(fc) - log(fm))/den; }
    return sigma;
  }

double neuromat_filter_lowpass_sgerf_compute_fsup(double fm, double sigma)
  { demand(isfinite(fm) && (fm > 0), "invalid {fm}");
    demand(isfinite(sigma) && (sigma > 0), "invalid {sigma}");
    return fm*exp(6.0*sigma);
  }

double neuromat_filter_lowpass_sgerf_eval(double f, double fm, double sigma)
  {
    demand(isfinite(fm) && (fm > 0), "invalid {fm}");
    demand(isfinite(sigma) && (sigma > 0), "invalid {sigma}");
    f = fabs(f);
    if (f == 0) { return 1.0; }
    double z = (log(f) - log(fm))/sigma;
    if (z >= 6.0) 
      { return 0; }
    else if (z <= -6.0)
      { return 1; }
    else
      { return 0.5*(1 - erf(z)); }
  }
