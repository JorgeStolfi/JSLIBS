#ifndef float_image_from_uint16_image_H
#define float_image_from_uint16_image_H

/* Conversion from multichannel {uint16_t}-valued images to float-valued images. */
/* Last edited on 2025-01-30 08:07:46 by stolfi */ 

#include <stdio.h>

#include <bool.h>
#include <jspnm.h>
#include <uint16_image.h>
#include <float_image.h>

float_image_t *float_image_from_uint16_image
  ( uint16_image_t *pim, /* PNM image to convert. */
    bool_t isMask,    /* TRUE for masks, FALSE for images. */
    double lo[],      /* Min float sample value. */
    double hi[],      /* Max float sample value. */
    bool_t yUp,       /* If TRUE, reverses the indexing of rows. */
    bool_t verbose    /* If TRUE, prints conversion diagnostics to {stderr}. */
  );
  /* Converts a PGM or PPM image {pim} into a {float_image_t} {fim}.
    Each sample {iv} in channel {k} of {pim} is converted to a float 
    with {sample_conv_floatize(iv, pim->maxval, isMask, lo[k], hi[k],...)}.
    
    If {lo} is NULL, it defaults to a vector of zeros.
    If {hi} is NULL, it defaults to a vector of ones.

    If {yUp} is TRUE, the row indices are reversed, so that row 0 of
    the resulting image {fim} is at the BOTTOM of the image; and, in
    general, row {i} of {fim} is row {ny-1-i} of {pim}. If {yUp} is
    FALSE, the row indices are preserved. NOTE: The {yUp} parameter was
    added on 2009-02-24.

    As a rough rule, the {isMask} flag should be FALSE for
    physically-based images (where the original sample values are
    expected to be smoothly distributed in {(lo_hi)}) and TRUE for
    mask images (where 0 and maxval must be mapped {lo} and {hi}
    exactly).  NOTE: The {isMask} parameter was added on 2010-08-14.

    If {verbose} is TRUE, the procedure prints statistics of the
    conversion to {stderr}. */

#endif
