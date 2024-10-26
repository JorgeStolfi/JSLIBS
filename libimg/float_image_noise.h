#ifndef float_image_noise_H
#define float_image_noise_H

/* Tools for generating images that are white, pink, or blue noise. */
/* Last edited on 2024-10-26 04:26:56 by stolfi */ 

#define _GNU_SOURCE_
#include <stdint.h>

#include <float_image.h>

float_image_t *float_image_noise
  ( int32_t NC,
    int32_t NX,
    int32_t NY,
    double fxFilter,
    double fyFilter,
    bool_t complement
  );
  /* Computes a noise image with given dimensions whose Hartley spectrum,
    in each channel, has the same amplitude at every non-zero frequency, 
    with random phase.
    
    The {fxFilter} and {fyFilter} parameters must be non-negative. If
    {fxFilter} is finite and positive, it specifies a Gaussian-like
    filter function {HX(fx) = exp(-(fx/fxFilter)^2/2)}, with foldover
    modulo {NX}. If {fxFilter} is zero ir {+INF}, the filter {HX(fx)} is
    1.0 (no filtering) for all {fx}.
    
    The {fyFilter} specifies an analogous filter {HY} for the vertical
    frequency {fy}. If {complement} is false, the amplitude of each
    component with integer frequency {(fx,fy)} is multiplied by by
    {HX(fx)*HY(fy)}. If {complement} is true, it is multiplied by
    {1-HX(fx)*HY(fy)}. Note that, in this case, if both {HX} and {HY}
    are 1.0 the result will be all zeros. */

#endif
