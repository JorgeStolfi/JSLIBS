/* See {neuromat_filter_bandpass_butterworth.h}. */
/* Last edited on 2024-01-05 18:16:13 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <jsmath.h>

#include <neuromat_filter_lowpass_butterworth.h>

#include <neuromat_filter_bandpass_butterworth.h>

void neuromat_filter_bandpass_butterworth_compute_parms
  ( double flo0,  double flo1,
    double fhi1,  double fhi0,
    double *fs_lo_P, double *fs_hi_P, int32_t *ord_P,
    double *fsup_P,
    bool_t verbose
  )
  {
    demand((!isnan(flo0)) && (!isnan(flo1)) && (!isnan(fhi1)) && (!isnan(fhi0)), "invalid NAN arguments");
    demand(flo1 < fhi1, "cannot have {flo1 >= fhi1}");
    demand(flo0 >= 0, "invalid {flo0}");
    double tiny = 1.0e-6;

    double fs_lo; int32_t ord_lo; /* Parameters of filter {W_lo}. */
    if (flo0 == 0)
      { demand(flo1 == 0, "{flo1} must be zero when {flo0} is zero");
        fs_lo = 0; ord_lo = 0;
      }
    else
      { demand((0 < flo0) && (flo0 < flo1), "invalid {flo0,flo1}"); 
        assert(flo1 < +INF); /* Since {flo1 < fhi1}. */
        fs_lo = flo1;
        ord_lo = neuromat_filter_lowpass_butterworth_compute_order(fs_lo, flo0, tiny);
      }
    if (verbose) 
      { fprintf(stderr, "  %s: flo0 = %16.12f flo1 = %16.12f", __FUNCTION__, flo0, flo1);
        fprintf(stderr, "  -->  fs_lo = %16.12f ord_lo = %d\n", fs_lo, ord_lo);
      }

    double fs_hi; int32_t ord_hi; /* Parameters of filter {W_hi}. */
    if (fhi0 == +INF)
      { demand(fhi1 == +INF, "{fhi1} must be {+INF} when {fhi0} is {+INF}");
        fs_hi = +INF; ord_hi = 0;
      }
    else
      { demand((fhi1 < fhi0) && (fhi0 < +INF), "invalid {fhi0,fhi1}"); 
        assert(fhi1 > 0); /* Since {flo1 < fhi1}. */
        fs_hi = fhi1;
        ord_hi = neuromat_filter_lowpass_butterworth_compute_order(fs_hi, fhi0, tiny);
      }
    if (verbose) 
      { fprintf(stderr, "  %s: fhi1 = %16.12f fhi0 = %16.12f", __FUNCTION__, fhi1, fhi0);
        fprintf(stderr, "  -->  fs_hi = %16.12f ord_hi = %d\n", fs_hi, ord_hi);
      }
    
    /* Take the max order. Note that it works even if either half-filter is suppressed. */
    int32_t ord = (ord_lo > ord_hi ? ord_lo : ord_hi);
      
    double fsup = fs_hi * exp(37.0/ord);
    if (verbose) { fprintf(stderr, "  %s: fsup = %16.12f\n", __FUNCTION__, fsup); }
  
    (*fs_lo_P) = fs_lo;
    (*fs_hi_P) = fs_hi;
    (*ord_P) = ord;
    (*fsup_P) = fsup;
  }    

double neuromat_filter_bandpass_butterworth_eval(double f, double fs_lo, double fs_hi, int32_t ord)
  {
    double W_lo;
    demand(isfinite(fs_lo) && (fs_lo >= 0), "invalid {fs_lo}");
    if (fs_lo == 0)
      { W_lo = 0; }
    else
      { W_lo = neuromat_filter_lowpass_butterworth_eval(f,fs_lo,ord); }
      
    double W_hi;
    demand(!isnan(fs_hi), "invalid {fs_hi}");
    if (fs_hi == +INF)
      { W_hi = 1; }
    else
      { W_hi = neuromat_filter_lowpass_butterworth_eval(f,fs_hi,ord); }
    
    return (1 - W_lo)*W_hi;
  }
