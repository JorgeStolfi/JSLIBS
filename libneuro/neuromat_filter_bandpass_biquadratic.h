#ifndef neuromat_filter_bandpass_biquadratic_H
#define neuromat_filter_bandpass_biquadratic_H

/* NeuroMat bandpass filter gain function with biquadratic shoulders in log freq. */
/* Last edited on 2024-01-05 17:35:38 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <bool.h>

double neuromat_filter_bandpass_biquadratic_eval(double f, double flo0, double flo1, double fhi1, double fhi0);
  /* Evaluates a bandpass filter function {W} at frequency {f}, whose
    shoulders are biquadratic sigmoids in log freq space, that are zero
    when {f <= flo0} or {f >= fhi0}, and 1 when {flo1 <= f <= fhi1}.
    
    Namely, the filter is {W(f) = (1-W_lo(f))*W_hi(f)} where
    {W_lo(f)=neuromat_filter_lowpass_biquadratic_eval(f,flo0,flo1)} and
    {W_hi(f)=neuromat_filter_lowpass_biquadratic_eval(f,fhi1,fhi0)}. The sign
    of {f} is ignored, so that {W} is even (symmetric about zero).
    
    If {flo0} is zero, {W_lo(f)} is always zero so {W} is the same as
    the low-pass filter {W_hi}. Symmetrically, if {fhi0} is {+INF},
    {W_hi(f)} is always 1 so that {W} is the same as the high-pass filter
    {1-W_lo}. */

#endif
