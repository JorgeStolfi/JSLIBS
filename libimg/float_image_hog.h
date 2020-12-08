#ifndef float_image_align_H
#define float_image_align_H

/* Histogram of gradient orientations for a float-valued image. */
/* Last edited on 2012-01-10 03:18:50 by stolfilocal */ 

#include <bool.h>
#include <float_image.h>

void float_image_hog_collect
  ( float_image_t *DX, 
    int cX,
    float_image_t *DY, 
    int cY,
    float_image_t *M, 
    double noise, 
    bool_t oriented,
    int nh,
    double h[]
  );
  /* Collects a histogram {h[0..nh-1]} of orientations of a gradient
    image, consisting of channel {cX} of image {DX} with channel {cY} of
    image {DY}. Assumes that each image is contaminated with
    Gaussian-like errors with deviation {noise}.
    
    The histogram consists of {nh} bins. Each pixel is added to the
    histogram, spread over one or more consecutive bins according to the
    direction uncertainty.  If {M} is not null, the total mass added for
    each pixel of the gradient is the value of the corresponding pixel 
    of {M}; otherwise the total mass is {1-exp((d/noise)^2/2)} where
    {d} is the modulus of the gradient. In particular, the pixel is ignored
    if {d} is zero.  The histogram is not normalized.
    
    Diametrally opposite gradient directions are considered distinct if
    {oriented} is true, identical if {orietned} is false. In the first
    case the histogram covers the full range {[0 _ 2*PI]}, in the second
    case the range {[0 _ PI]}. The first bin is centered at angle zero
    (direction {(1,0)}) and move on counterclockwise. */

#endif
