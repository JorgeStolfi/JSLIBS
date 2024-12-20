#ifndef float_array_to_uint16_image_H
#define float_array_to_uint16_image_H

/* Conversion between {uint16_image_t} images and multi-dimensional float arrays. */
/* Last edited on 2024-12-05 10:31:55 by stolfi */ 

#include <stdint.h>
#include <stdio.h>

#include <bool.h>
#include <uint16_image.h>

#include <float_array.h>

float_array_t *float_array_from_uint16_image
  ( uint16_image_t *img, 
    bool_t isMask,
    double lo[], 
    double hi[], 
    bool_t verbose
  );
  /* Converts a PGM or PPM image into a {float_array_t}.
    Each sample {iv} in channel {k} is converted to a float 
    with {sample_conv_floatize(iv,img->maxval,isMask,lo[k],hi[k],...)}
    (see {sample_conv.h}).

    If the vector {lo} is NULL, it defaults to a vector of zeros.
    If the vector {hi} is NULL, it defaults to a vector of ones.
    
    If {img} is grayscale, the channel index {ix[0]} has only one
    legal value, namely 0; and each sample is the intensity (Y) of the
    pixel. If {img} is a color image, then {ix[0] = 0,1,2} select the
    red, green, and blue channel, respectively.

    If {verbose} is TRUE, the procedure prints statistics of the
    conversion to {stderr}. */

uint16_image_t *float_array_to_uint16_image
  ( float_array_t *A,
    bool_t isMask,
    int32_t chns,
    int32_t ch[],
    double lo[], 
    double hi[], 
    bool_t yrev,
    uint16_t maxval, 
    bool_t verbose
  );
  /* Converts a {float_array_t} {A} (with float samples) to
    a {uint16_image_t} {img} (with integer samples).  

    The image {img} will have {chns} channels, numbered {0..chns-1} The
    samples of channel {k} of {img} will be taken from the subarray of
    {A} with {ix[0] = ch[k]}, for each {k} in {0..chns-1}. If the index
    {ch[k]} is invalid (not in {0..A.sz[0]-1}), the
    samples in channel {k} of {img} are all set to zero. If the {ch}
    argument is NULL, it defaults to the identity vector
    {(0,1,2,...chns-1)}.
    
    The image will have {img.cols=A.sz[1]} and {img.rows=A.sz[2]}.
    The numbering of columns will be preserved.
    The numbering of the rows will be preserved if
    {yrev} is false, or reversed if {yrev} is true.

    The remaining indices of {A} must be trivial trivial;
    that is, {A.sz[k]=1} for {k} in {0..float_array_MAX_AXES-1}.
    
    Each element {fv} of {A} destined to be a sample in channel {k} of
    {img} is converted to an integer sample in {0..maxval} by the call
    {sample_conv_quantize(fv,maxval,isMask,lo[c],hi[c],...)} (see
    {sample_conv.h}). If {lo} is NULL, it defaults to a vector of zeros.
    If {hi} is NULL, it defaults to a vector of ones.

    If {verbose} is TRUE, the procedure prints statistics of the
    conversion to {stderr}. */

#endif
