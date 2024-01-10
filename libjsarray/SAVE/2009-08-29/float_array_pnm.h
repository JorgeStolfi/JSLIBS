#ifndef float_array_pnm_H
#define float_array_pnm_H

/* Conversion between PPM images and multi-dimensional float images. */
/* Last edited on 2007-10-11 01:06:54 by stolfi */ 

#include <stdio.h>

#include <bool.h>
#include <jspnm.h>
#include <jspnm_image.h>

#include <float_array.h>

/* CONVERSION BETWEEN PPM/PBM/PGM IMAGES AND FLOAT ARRAYS

  The procedures in this interface perform conversion between
  a {float_array_t} {A} (with float samples) and a {pnm_image_t} {img}
  (with integer samples).  

  The indices of the {float_array_t} will have the following meaning:

    | {ix[0]} = color channel.
    | {ix[1]} = pixel column (leftmost=0, rightmost = {width - 1}).
    | {ix[2]} = pixel row (bottom = 0, top = {height - 1}).

  The remaining indices {ix[3..5]} are trivial (always 0).
  
  For the meaning of the color channel index {ix[0]},see each routine.

  Note that the values of the row index {ix[2]} are reversed with
  respect to those of {img}. Namely, row 0 of {img} (the topmost one)
  is row {height-1} of {A}, and vice-versa. */

float_array_t float_array_from_pnm_image
  ( pnm_image_t *img, 
    double lo[], 
    double hi[], 
    bool_t verbose
  );
  /* Converts a PGM or PPM image into a {float_array_t}.
    Each sample {iv} in channel {k} is converted to a float 
    with {sample_conv_floatize(iv,img->maxval,lo[k],hi[k],...)}
    (see {sample_conv.h}).

    If the vector {lo} is NULL, it defaults to a vector of zeros.
    If the vector {hi} is NULL, it defaults to a vector of ones.
    
    If {img} is grayscale, the channel index {ix[0]} has only one
    legal value, namely 0; and each sample is the intensity (Y) of the
    pixel. If {img} is a color image, then {ix[0] = 0,1,2} select the
    red, green, and blue channel, respectively.

    If {verbose} is TRUE, the procedure prints statistics of the
    conversion to {stderr}. */

pnm_image_t *float_array_to_pnm_image
  ( float_array_t *A,
    int chns,
    int ch[],
    double lo[], 
    double hi[], 
    pnm_sample_t maxval, 
    bool_t verbose
  );
  /* Converts a {float_array_t} {A} to a PGM/PPM image {pim}.
  
    The number of channels of {pim} with be {chns}, which is usually 
    1 for grayscale (PGM) and 3 for color (PPM). 
    
    The samples of channel {k} of {pim} will be taken from the
    subarray of {A} with {ix[0] = ch[k]}, for each {k} in {0..chns-1}.
    If {ch[k]} is invalid (negative, or {>= A.sz[0]}), the samples in
    channel {k} of {img} are all set to zero. If {ch} is NULL, it
    defaults to the identity vector {(0,1,2,...chns-1)}.
    
    Each sample {fv} destined to channel {k} of {img} is converted to
    an integer sample in {0..maxval} by the call
    {sample_conv_quantize(fv,maxval,lo[c],hi[c],...)} (see
    {sample_conv.h}). If {lo} is NULL, it defaults to a vector of
    zeros. If {hi} is NULL, it defaults to a vector of ones.
    
    Note that the row indices are reversed, so that row 0 of {img}
    will be row {img->rows-1} of the result, and vice-versa.

    If {verbose} is TRUE, the procedure prints statistics of the
    conversion to {stderr}. */

#endif
