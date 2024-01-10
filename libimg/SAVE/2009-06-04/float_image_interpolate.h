#ifndef float_image_interpolate_H
#define float_image_interpolate_H

/* Bilinear (C0) and Bicubic (C1) interpolation of floating-point images. */
/* Last edited on 2009-06-04 10:49:33 by stolfi */ 

#include <bool.h>
#include <indexing.h>
#include <float_image.h>

/* 
  PIXEL INDICES AND COORDINATES
  
  The procedures in this module assume that the pixel on column {ix} and 
  row {iy} is a unit square centered at point {(ix+0.5,iy+0.5)}.  
  
  PIXELS OUTSIDE THE IMAGE DOMAIN
  
  Pixel indices outside the valid ranges {0..NX-1}, {0..NY-1} are
  reduced according to {ix_reduce(i, N, red)}. If a reduced pixel
  index is {-1}, the functions assume that its sample values are all 0.
  
  INTERPOLATION METHODS
  
  If the {order} parameter is -1, the interpolation is piecewise constant
  and discontinuous: the 
  result is the value of the pixel that contains the given point {(x,y)}.
  
  If {order} is 0, the interpolation is piecewise affine and C0. The
  result is a weighted average of the 2×2 block of pixels whose
  centers surround the point {(x,y)}. It reconstructs bilinear
  functions.
  
  If {order} is 1, the interpolation is bicubic and C1. The result is
  a weighted average of the 4×4 block of pixels such that the point
  {(x,y)} is surrounded by the centers of the innermost 2×2 block of
  pixels. It reconstructs biquadratic functions.  */

double float_image_interpolate_sample(float_image_t *A, int c, double x, double y, int order, ix_reduction_t red);
  /* Returns the value of channel {c} of image {A} at the point
    {(x,y)}, computed by interpolation from nearby samples. */

void float_image_interpolate_pixel(float_image_t *A, double x, double y, int order, ix_reduction_t red, double v[]);
  /* For {c=0..NC-1}, stores into {v[c]} the value of chanel {c} of
    image {A} at the point {(x,y)}, computed as in
    {float_image_interpolate_C1_sample}; where {NC} is the number of
    channels of {A}. */

#endif
