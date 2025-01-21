#ifndef float_image_shrink_by_one_H
#define float_image_shrink_by_one_H

/* Shrink an image by one row and one column. */
/* Last edited on 2025-01-18 17:35:33 by stolfi */ 

#include <float_image.h>

float_image_t *float_image_shrink_by_one(float_image_t *A, int32_t cW);
  /* Shrinks the image {A} by one col and one row. In the result {R},
    each sample {R[c,x,y} is the average of the four samples
    {R[c,x+dx,y+dy]} where {dx,dy} are 0 or 1.
    
    If {cW} is the index of a channel of {A}, sample {A[cW,x,y]} is
    assumed to be the weight sample {A[c,x,y]}, for each {c != cW}, for
    the purpose of averaging. Otherwise the weight is assumed to be 1
    for all pixels of {A}. Sample {R[cW,x,y]} is then set to the
    (unweighted) harmonic mean of the weights of the pixels used to
    compute {R[c,x,y]}. In particular, if any of those weights is zero,
    {R[cW,x,y]} will be zero. */

#endif
