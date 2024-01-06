/* See {neuromat_filter_bandpass_biquadratic.h}. */
/* Last edited on 2024-01-05 17:35:50 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <jsmath.h>

#include <neuromat_filter_lowpass_biquadratic.h>

#include <neuromat_filter_bandpass_biquadratic.h>

double neuromat_filter_bandpass_biquadratic_eval(double f, double flo0, double flo1, double fhi1, double fhi0)
  {
    double W_lo = (flo0 == 0 ? 0 : neuromat_filter_lowpass_biquadratic_eval(f,flo0,flo1));
    double W_hi = (fhi0 == +INF ? 1 : neuromat_filter_lowpass_biquadratic_eval(f,fhi1,fhi0));
    return (1 - W_lo)*W_hi;
  }
