#ifndef neuromat_filter_lowpass_biquadratic_H
#define neuromat_filter_lowpass_biquadratic_H

/* NeuroMat basic lowpass filter transfer function with biquadratic shoulder in log frequency space. */
/* Last edited on 2024-01-05 17:33:48 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
    
double neuromat_filter_lowpass_biquadratic_eval(double f, double fa, double fb);
  /* A filter that rolls off from gain 1 at {f==fa} to gain zero at {f==fb}.
    The filter is a sigmoid-like biquadratic in log frequency space.
    Requires {0 <= fa < fb}.  
    
    If {f} is zero, returns 1.  Else, if {fa} is zero, returns zero.
    If {f} is negative, returns the same value as for {-f}. */ 

#endif
