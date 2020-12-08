#ifndef float_image_write_pnm_H
#define float_image_write_pnm_H

/* Writing float images to PNM (PPM/PGM/PBM) image files. */
/* Last edited on 2017-06-25 00:14:15 by stolfilocal */

#define _GNU_SOURCE
#include <bool.h>
#include <float_image.h>

void float_image_write_pnm_named
  ( char *fname,        /* PPM/PGM/PBM file name (with extension). */            
    float_image_t *fim, /* Image to write. */
    bool_t isMask,      /* TRUE for masks, FALSE for images. */
    double gamma,       /* Gamma to use in encoding (1 = linear encoding). */    
    double bias,        /* Offset to use in encoding. */                         
    bool_t yup,         /* If TRUE, reverses the indexing of rows. */ 
    bool_t warn,        /* If TRUE, prints "writing {fname}..." to {stderr}. */
    bool_t verbose      /* If TRUE, prints conversion diagnostics to {stderr}. */
  );
  /* Writes the image as a PGM or PPM image file, depending on the
    number of channels. Clips the given image samples to [0_1] and
    quantizes them with {maxval=uint16_image_MAX_SAMPLE}. Applies the gamma
    encoding corresponding to the decoding exponent {gamma} and offset {bias}.
   
    See {sample_conv_floatize} for the meaning of the {isMask} parameter.
    NOTE: The{isMask} parameter was added on 2010-08-14.

    If {yup} is TRUE, reverses the indices of rows, so that row 0 will
    be at the BOTTOM of the image. NOTE: This parameter was
    added on 2009-02-24.
    
    If {warn} is true, prints "writing {fname}..." to {stderr}. If
    {verbose} is true, prints also the sample conversion statistics.
    NOTE: These parameters were added on 2009-02-24. */

void float_image_write_pnm_named_list
  ( int n,                /* Number of images to write. */
    char *fname[],        /* PPM/PGM/PBM file names (with extensions). */            
    float_image_t *fim[], /* Images to write. */
    bool_t isMask,        /* TRUE for masks, FALSE for images. */
    double gamma,         /* Gamma to use in encoding (1 = linear encoding). */    
    double bias,          /* Offset to use in encoding. */                         
    bool_t yup,           /* If TRUE, reverses the indexing of rows. */ 
    bool_t warn,          /* If TRUE, prints "writing {fname}..." to {stderr}. */
    bool_t verbose        /* If TRUE, prints conversion diagnostics to {stderr}. */
  );
  /* Calls {float_image_write_pnm_named_pnm_named(fname[i],fim[i],isMask,gamma,bias,yup,warn,verbose)}
    for {i} in {0..n-1}. 
    
    NOTE: order of {n,fname} changed and {yp,warn,verbose} added on
    2009-02-24. */

#endif
