/* See {neuromat_filter_bandpass_sgerf.h}. */
/* Last edited on 2024-01-06 08:17:56 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <jsmath.h>

#include <neuromat_filter_lowpass_sgerf.h>

#include <neuromat_filter_bandpass_sgerf.h>

void neuromat_filter_bandpass_sgerf_compute_parms
  ( double flo0,  double flo1,
    double fhi1,  double fhi0,
    double *fm_lo_P, double *sigma_lo_P,
    double *fm_hi_P, double *sigma_hi_P,
    double *fsup_P,
    bool_t verbose
  )
  { demand((!isnan(flo0)) && (!isnan(flo1)) && (!isnan(fhi1)) && (!isnan(fhi0)), "invalid NAN arguments");
    demand(flo1 < fhi1, "cannot have {flo1 >= fhi1}");
    demand(flo0 >= 0, "invalid {flo0}");
    double tiny = 1.0e-5;

    double fm_lo, sigma_lo; /* Parameters of filter {W_lo}. */
    if (flo0 == 0)
      { demand(flo1 == 0, "{flo1} must be zero when {flo0} is zero");
        fm_lo = 0; sigma_lo = 0;
      }
    else
      { demand((0 < flo0) && (flo0 < flo1), "invalid {flo0,flo1}"); 
        assert(flo1 < +INF); /* Since {flo1 < fhi1}. */
        fm_lo = sqrt(flo0*flo1);
        sigma_lo = neuromat_filter_lowpass_sgerf_compute_sigma(fm_lo, flo0, tiny);
      }
    if (verbose) 
      { fprintf(stderr, "  %s: flo0 = %16.12f flo1 = %16.12f", __FUNCTION__, flo0, flo1);
        fprintf(stderr, "  -->  fm_lo = %16.12f sigma_lo = %16.12f\n", fm_lo, sigma_lo);
      }

    double fm_hi, sigma_hi; /* Parameters of filter {W_hi}. */
    if (fhi0 == +INF)
      { demand(fhi1 == +INF, "{fhi1} must be {+INF} when {fhi0} is {+INF}");
        fm_hi = +INF; sigma_hi = 0;
      }
    else
      { demand((fhi1 < fhi0) && (fhi0 < +INF), "invalid {fhi0,fhi1}"); 
        assert(fhi1 > 0); /* Since {flo1 < fhi1}. */
        fm_hi = sqrt(fhi0*fhi1);
        sigma_hi = neuromat_filter_lowpass_sgerf_compute_sigma(fm_hi, fhi0, tiny);
      }
    if (verbose) 
      { fprintf(stderr, "  %s: fhi1 = %16.12f fhi0 = %16.12f", __FUNCTION__, fhi1, fhi0);
        fprintf(stderr, "  -->  fm_hi = %16.12f sigma_hi = %16.12f\n", fm_hi, sigma_hi);
      }
      
    double fsup = exp(log(fm_hi) + 6*sigma_hi);
    if (verbose) { fprintf(stderr, "  %s: fsup = %16.12f\n", __FUNCTION__, fsup); }
        
    (*fm_lo_P) = fm_lo;
    (*sigma_lo_P) = sigma_lo;
    (*fm_hi_P) = fm_hi;
    (*sigma_hi_P) = sigma_hi;
    (*fsup_P) = fsup;
  }    

double neuromat_filter_bandpass_sgerf_eval(double f, double fm_lo, double sigma_lo, double fm_hi, double sigma_hi)
  {
    double W_lo;
    demand(isfinite(fm_lo) && (fm_lo >= 0), "invalid {fm_lo}");
    if (fm_lo == 0)
      { W_lo = 0; }
    else
      { W_lo = neuromat_filter_lowpass_sgerf_eval(f,fm_lo,sigma_lo); }
      
    double W_hi;
    demand(!isnan(fm_hi), "invalid {fm_hi}");
    if (fm_hi == +INF)
      { W_hi = 1; }
    else
      { W_hi = neuromat_filter_lowpass_sgerf_eval(f,fm_hi,sigma_hi); }
    
    return (1 - W_lo)*W_hi;
  }
