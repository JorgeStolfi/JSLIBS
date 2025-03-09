#ifndef float_image_expand_by_one_H
#define float_image_expand_by_one_H

/* Expand an image by one row and one column. */
/* Last edited on 2025-02-22 16:58:24 by stolfi */ 

#include <float_image.h>

float_image_t *float_image_expand_by_one(float_image_t *A, int32_t cW);
  /* Expands the image {A} by one col and one row. In the resulting
    image {E}, each pixel {E[c,x,y]} is the average of the pixels
    {E[c,x-dx,y-dy]} where {dx,dy} are 0 or 1, if those pixels exist.
    
    If {cW} is the valid index of a channel of {A}, sample {A[cW,x,y]} is
    assumed to be the weight sample {A[c,x,y]}, for each {c != cW}, for
    the purpose of averaging. Otherwise the weight is assumed to be 1
    for all pixels of {A}. Sample {E[cW,x,y]} is then set to the
    (unweighted) harmonic mean of the weights of the pixels used to
    compute {E[c,x,y]}.  In particular, if any of those weights is zero,
    {E[cW,x,y]} will be zero. */

#endif
