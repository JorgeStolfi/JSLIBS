#ifndef float_image_io_pnm_H
#define float_image_io_pnm_H

/* Reading/writing float images from/to PPM/PGM/PBM files. */
/* Last edited on 2017-06-20 20:48:41 by stolfilocal */

#include "float_image.h"

#define VIEW_GAMMA 2.200 
#define VIEW_BIAS  0.010
  /* Gamma and bias for viewing images on IBM-PC-like platforms. */

/* SINGLE IMAGE I/O */

float_image_t *float_image_from_to_uint16_image_read
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

void float_image_from_to_uint16_image_write
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
    quantizes them with {PNM_FILE_MAX_MAXVAL}. Applies the gamma
    encoding with decoding exponent {gamma} and offset {bias}.
   
    See {sample_conv_floatize} for the meaning of the {isMask} parameter.
    NOTE: The{isMask} parameter was added on 2010-08-14.

    If {yup} is TRUE, reverses the indices of rows, so that row 0 will
    be at the BOTTOM of the image. NOTE: This parameter was
    added on 2009-02-24.
    
    If {warn} is true, prints "writing {fname}..." to {stderr}. If
    {verbose} is true, prints also the sample conversion statistics.
    NOTE: These parameters were added on 2009-02-24. */

/* IMAGE LIST I/O */

float_image_t **float_image_from_to_uint16_image_list_read
  ( int n,          /* Number of images to read. */
    char *fname[],  /* PPM/PGM/PBM file names (with extensions). */
    bool_t isMask,  /* TRUE for masks, FALSE for images. */
    double gamma,   /* Gamma to use in decoding (1 = linear decoding). */
    double bias,    /* Offset to use in decoding. */
    bool_t yup,     /* If TRUE, reverses the indexing of rows. */
    bool_t warn,    /* If TRUE, prints "reading {fname}..." to {stderr}. */
    bool_t verbose  /* If TRUE, prints conversion diagnostics to {stderr}. */
  );
  /* Calls {float_image_from_to_uint16_image_read(fname[i],isMask,gamma,bias,yup,warn,verbose)}
     for {i} in {0..n-1} and returns an array of pointers to the resulting
     float images. 

     NOTE: order of {n,fname} changed and {yp,warn,verbose} added on
     2009-02-24. */

void float_image_from_to_uint16_image_list_write
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
  /* Calls {float_image_from_to_uint16_image_write(fname[i],fim[i],isMask,gamma,bias,yup,warn,verbose)}
    for {i} in {0..n-1}. 
    
    NOTE: order of {n,fname} changed and {yp,warn,verbose} added on
    2009-02-24. */

#endif
