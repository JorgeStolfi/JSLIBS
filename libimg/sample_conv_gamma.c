/* See {sample_conv_gamma.h}. */
/* Last edited on 2024-12-18 23:04:40 by stolfi */

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>

#include <affirm.h>
#include <bool.h>
#include <jsmath.h>
#include <sample_conv.h>

#include <sample_conv_gamma.h>

float sample_conv_gamma(float z, double expo, double bias)
  { demand(bias >= 0, "negative bias not implemented yet");
    if(isnan(z) || isinf(z)) { return z; }
    if (expo == 1) { return z; }
    float a = fabsf(z);
    if ((a == 0) || (a == 1)) { return z; }
    if (bias == 0)
      { a = (float)pow(a, expo); }
    else
      { double sg = sqrt(expo);
        double c = pow(bias, 1/sg);
        double d = pow(bias, sg);
        double u = a*(1 - c) + c;
        double v = pow(u, expo);
        double w = (v - d)/(1 - d);
        a = (float)w;
        if (a < 0) { a = 0; }
      }
    return (z < 0 ? -a : a);
  }
