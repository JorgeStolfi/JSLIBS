/* See {neuromat_filter_bandpass_gauss.h}. */
/* Last edited on 2024-01-05 18:15:41 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <jsmath.h>

#include <neuromat_filter_lowpass_gauss.h>

#include <neuromat_filter_bandpass_gauss.h>
 
void neuromat_filter_bandpass_gauss_compute_parms
  ( double flo0,  double flo1,
    double fhi1,  double fhi0,
    double *sigma_lo_P, double *sigma_hi_P, double *mag_P,
    double *fsup_P,
    bool_t verbose
  )
  {
    demand((!isnan(flo0)) && (!isnan(flo1)) && (!isnan(fhi1)) && (!isnan(fhi0)), "invalid NAN arguments");
    demand(flo0 >= 0, "invalid {flo0}");
    demand(flo1 < fhi1, "cannot have {flo1 >= fhi1}");
    double tiny = 1.0e-6;

    double sigma_hi = NAN;
    if (fhi0 == +INF)
      { demand(fhi1 == +INF, "{fhi1} must be {+INF} when {fhi0} is {+INF}");
        sigma_hi = +INF;
      }
    else
      { demand(fhi1 < fhi0, "invalid {fhi1, fhi0}");
        assert(fhi1 > 0); /* Since {flo1 < fhi1}. */
        sigma_hi = neuromat_filter_lowpass_gauss_compute_sigma(fhi0, tiny);
      }
    if (verbose) 
      { fprintf(stderr, "  %s: fhi1 = %16.12f fhi0 = %16.12f", __FUNCTION__, fhi1, fhi0);
        fprintf(stderr, "  -->  sigma_hi = %16.12f\n", sigma_hi);
      }

    double sigma_lo = NAN;
    if (flo0 == 0)
      { demand(flo1 == 0, "{flo1} must be zero when {flo0} is zero");
        sigma_lo = 0;
      }
    else
      { demand(flo0 < flo1, "invalid {fhi1, fhi0}");
        assert(flo1 < +INF); /* Since {flo1 < fhi1}. */
        sigma_lo = neuromat_filter_lowpass_gauss_compute_sigma(fhi0, 0.5);
        sigma_lo = fmin(sigma_lo, 0.5*sigma_hi);
      }
    if (verbose) 
      { fprintf(stderr, "  %s: flo0 = %16.12f flo1 = %16.12f", __FUNCTION__, flo0, flo1);
        fprintf(stderr, "  -->  sigma_lo = %16.12f\n", sigma_lo);
      }

    double fmid = (flo1 + fhi1)/2;
    double gmid = neuromat_filter_bandpass_gauss_eval(fmid, sigma_lo, sigma_hi, 1.0);
    double mag = 1.0/gmid;
    if (verbose) 
      { fprintf(stderr, "  %s: fmid = %16.12f gmid = %16.12f", __FUNCTION__, fmid, gmid);
        fprintf(stderr, "  -->  mag = %16.12f\n", mag);
      }
      
    double fsup = 9*sigma_hi;
    if (verbose) { fprintf(stderr, "  %s: fsup = %16.12f\n", __FUNCTION__, fsup); }
    (*fsup_P) = fsup;
  
    (*sigma_lo_P) = sigma_lo;
    (*sigma_hi_P) = sigma_hi;
    (*mag_P) = mag;
    (*fsup_P) = fsup;
  }    

double neuromat_filter_bandpass_gauss_eval(double f, double sigma_lo, double sigma_hi, double mag)
  {
    double W_lo;
    demand(isfinite(sigma_lo) && (sigma_lo >= 0), "invalid {sigma_lo}");
    if (sigma_lo == 0)
      { W_lo = 0; }
    else
      { W_lo = neuromat_filter_lowpass_gauss_eval(f,sigma_lo); }
      
    double W_hi;
    demand(!isnan(sigma_hi), "invalid {sigma_hi}");
    if (sigma_hi == +INF)
      { W_hi = 1; }
    else
      { W_hi = neuromat_filter_lowpass_gauss_eval(f,sigma_hi); }
      
    /* fprintf(stderr, "    %s: sigma_lo = %16.12f W_lo = %16.12f\n", __FUNCTION__, sigma_lo, W_lo);  */
    /* fprintf(stderr, "    %s: sigma_hi = %16.12f W_hi = %16.12f\n", __FUNCTION__, sigma_hi, W_hi);  */
   
    return mag*(1 - W_lo)*W_hi;
  }
