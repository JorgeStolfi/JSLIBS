#ifndef float_array_from_uint16_image_H
#define float_array_from_uint16_image_H

/* Conversion from {uint16_image_t} images to multi-dimensional float arrays. */
/* Last edited on 2024-12-05 10:31:50 by stolfi */ 

#include <stdio.h>

#include <bool.h>
#include <uint16_image.h>

#include <float_array.h>

float_array_t *float_array_from_uint16_image
  ( uint16_image_t *img, 
    bool_t isMask,
    double lo[], 
    double hi[], 
    bool_t yrev,
    bool_t verbose
  );
  /* Converts a {uint16_image_t} {img} (with integer samples) to a
    {float_array_t} {A} (with float samples) and .

    Indices {ix[0], ix[1], and ix[2]} of the {float_array_t} will correstpon to the
    indices of channel, pixel column, and pixel row of the image {img},
    and will have the same ranges -- that is, {A.sz[0]=img.chns}, {A.sz[1]=img.cols},
    and {A.sz[2]=img.rows}. . 
    
    The numbering of channels and columns will be preserved.
    The numbering of the rows will be preserved if
    {yrev} is false, or reversed if {yrev} is true.

    The remaining indices of {A} will be trivial, namely 
    {A.sz[k]=1} for {k} in {0..float_array_MAX_AXES-1}.
    
    Each sample {iv} in channel {k} is converted to a float 
    with {sample_conv_floatize(iv,img->maxval,isMask,lo[k],hi[k],...)}
    (see {sample_conv.h}).

    If the vector {lo} is NULL, it defaults to a vector of zeros.
    If the vector {hi} is NULL, it defaults to a vector of ones.

    If {verbose} is TRUE, the procedure prints statistics of the
    conversion to {stderr}. */

#endif
