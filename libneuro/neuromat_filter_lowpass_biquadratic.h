#ifndef neuromat_filter_lowpass_biquadratic_H
#define neuromat_filter_lowpass_biquadratic_H

/* NeuroMat basic lowpass filter transfer function with biquadratic shoulder in log frequency space. */
/* Last edited on 2023-12-15 13:25:39 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
    
double neuromat_filter_lowpass_biquadratic(double f, double fa, double fb);
  /* A filter that rolls off from gain 1 at {f==fa} to gain zero at {f==fb}.
    The filter is a sigmoid-like biquadratic in log frequency space.
    Requires {0 <= fa < fb}. If {fa} is zero, returns zero. */ 

#endif
