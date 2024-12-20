#ifndef float_image_read_pnm_H
#define float_image_read_pnm_H

/* Reading float images from PNM (PPM/PGM/PBM) image files. */
/* Last edited on 2024-12-05 10:29:54 by stolfi */

#include <bool.h>
#include <float_image.h>

float_image_t *float_image_read_pnm_named
  ( char *fname,    /* PPM/PGM/PBM file name (with extension). */
    bool_t isMask,  /* TRUE for masks, FALSE for images. */
    double gamma,   /* Gamma to use in decoding (1 = linear decoding). */
    double bias,    /* Offset to use in decoding. */
    bool_t yup,     /* If TRUE, reverses the indexing of rows. */
    bool_t warn,    /* If TRUE, prints "reading {fname}..." to {stderr}. */
    bool_t verbose  /* If TRUE, prints conversion diagnostics to {stderr}. */
  );
  /* Reads the PGM or PPM image file with the given {fname}, and
    converts it to a float image with samples in [0_1]. Applies the
    gamma decoding with exponent {gamma} and offset {bias}.

    See {sample_conv_floatize} for the meaning of the {isMask} parameter.
    NOTE: The{isMask} parameter was added on 2010-08-14.

    If {yup} is TRUE, reverses the indices of rows, so that row 0 will
    be at the BOTTOM of the image. NOTE: This parameter was
    added on 2009-02-24.

    If {warn} is TRUE, prints "reading {fname}..." to {stderr}. If
    {verbose} is true, prints also the sample conversion statistics.
    NOTE: These parameters were added on 2009-02-24. */

float_image_t **float_image_read_pnm_named_list
  ( int32_t n,          /* Number of images to read. */
    char *fname[],  /* PPM/PGM/PBM file names (with extensions). */
    bool_t isMask,  /* TRUE for masks, FALSE for images. */
    double gamma,   /* Gamma to use in decoding (1 = linear decoding). */
    double bias,    /* Offset to use in decoding. */
    bool_t yup,     /* If TRUE, reverses the indexing of rows. */
    bool_t warn,    /* If TRUE, prints "reading {fname}..." to {stderr}. */
    bool_t verbose  /* If TRUE, prints conversion diagnostics to {stderr}. */
  );
  /* Calls {float_image_read_pnm_named_pnm_named(fname[i],isMask,gamma,bias,yup,warn,verbose)}
     for {i} in {0..n-1} and returns an array of pointers to the resulting
     float images. 

     NOTE: order of {n,fname} changed and {yp,warn,verbose} added on
     2009-02-24. */

#endif
