/* See {neuromat_filter_lowpass_biquadratic.h}. */
/* Last edited on 2023-12-15 13:26:52 by stolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <jsmath.h>

#include <neuromat_filter_lowpass_biquadratic.h>

double neuromat_filter_lowpass_biquadratic(double f, double fa, double fb)
  { demand((0 <= fa) && (fa < fb), "invalid {fa, fb}");
    if (fa == 0) { return 0.0; }
    if (fabs(f) <= fa) { return 1.0; }
    if (fabs(f) >= fb) { return 0.0; }
    double z = (log(fabs(f)) - log(fa))/(log(fb) - log(fa));
    double y = 1 - z*z;
    double g = y*y;
    return g;
  } 
