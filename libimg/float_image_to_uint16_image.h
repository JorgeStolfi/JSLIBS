#ifndef float_image_to_uint16_image_H
#define float_image_to_uint16_image_H

/* Conversion from multichannel float-valued images to {uint16_t}-valued images. */
/* Last edited on 2017-06-22 17:41:48 by stolfilocal */ 

#include <stdio.h>

#include <bool.h>
#include <jspnm.h>
#include <uint16_image.h>
#include <float_image.h>

uint16_image_t *float_image_to_uint16_image
  ( float_image_t *fim,  /* Float image to convert. */
    bool_t isMask,       /* TRUE for masks, FALSE for images. */
    int chns,            /* Number of channels of output image. */
    double lo[],         /* Nominal min float sample for each chosen channel. */
    double hi[],         /* Nominal max float sample for each chosen channel. */
    int ch[],            /* Indices of channels of {fim} to convert. */
    uint16_t maxval, /* Max integer sample value in result image. */
    bool_t yup,          /* If TRUE, reverses the indexing of rows. */
    bool_t verbose       /* If TRUE, prints conversion diagnostics to {stderr}. */
  );
  /* Converts a {float_image_t} {fim} to a PGM/PPM image {pim}. The
    number of channels of {pim} will be {chns}, and the samples of
    channel {k} of {pim} will be taken from channel {ch[k]} of {fim}.
    If a channel index {ch[k]} is invalid, its samples are assumed to
    be all zeros.  If {ch} is NULL, it defaults to the identity vector
    {(0,1,2,...chns-1)}.
    
    Each sample {fv} of channel {c = ch[k]} of {fim} is converted to
    an integer sample in {0..maxval} by {float_image_quantize(fv,
    maxval, isMask, lo[k], hi[k], ...)}. If {lo} is NULL, it defaults to a
    vector of zeros. If {hi} is NULL, it defaults to a vector of ones.
    
    If {yup} is TRUE, the row indices are reversed, so that 
    row 0 of the {float_image_t} is at the BOTTOM of the image;
    and, in general, row {i} of the result is row {ny-1-i} of {pim}. 
    If {yup} is FALSE, the row indices are preserved. 
    NOTE: The {yup} parameter was added on 2009-02-24.

    As a rough rule, the {isMask} flag should be FALSE for
    physically-based images (where the original sample values are
    expected to be smoothly distributed in {(lo_hi)}) and TRUE for
    mask images (where 0 and maxval must be mapped {lo} and {hi}
    exactly). NOTE: The{isMask} parameter was added on 2010-08-14.

    If {verbose} is TRUE, the procedure prints statistics of the
    conversion to {stderr}.
    
    NOTE: the indexing of {lo} and {hi} changed from {c} 
    to {k} on 2008-11-10.  The order of the parameters
    {lo,hi,ch} was changed to make old uses stand out. */

#endif
