/* See {float_image_write_pnm.h}. */
/* Last edited on 2024-12-04 23:32:49 by stolfi */

#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <float_image.h>
#include <uint16_image.h>
#include <uint16_image_write_pnm.h>
#include <float_image_to_uint16_image.h>

#include <float_image_write_pnm.h>

void float_image_write_pnm_named
  ( char *fname,        /* PPM/PGM/PBM file name (with extension). */            
    float_image_t *fim, /* Image to write. */
    bool_t isMask,      /* TRUE for masks, FALSE for images. */
    double gamma,       /* Gamma to use in encoding (1 = linear encoding). */    
    double bias,        /* Offset to use in encoding. */                         
    bool_t yup,         /* If TRUE, reverses the indexing of rows. */            
    bool_t warn,        /* If TRUE, prints "writing {fname}..." to {stderr}. */
    bool_t verbose      /* If TRUE, prints conversion diagnostics to {stderr}. */
  )
  { int32_t nc = (int32_t)fim->sz[0];
    assert((nc == 1) || (nc == 3));
    float_image_t * gim = float_image_copy(fim);
    int32_t c;
    for (c = 0; c < fim->sz[0]; c++)
      { float_image_apply_gamma(gim, c, 1.0/gamma, bias); }
    uint16_t maxval = uint16_image_MAX_SAMPLE;
    uint16_image_t *pim = 
      float_image_to_uint16_image(gim, isMask, nc, NULL, NULL, NULL, maxval, yup, verbose);
    float_image_free(gim);
    uint16_image_write_pnm_named(fname, pim, FALSE, warn);
    uint16_image_free(pim);
  }

void float_image_write_pnm_named_list
  ( int32_t n,                /* Number of images to write. */
    char *fname[],        /* PPM/PGM/PBM file names (with extensions). */            
    float_image_t *fim[], /* Images to write. */
    bool_t isMask,        /* TRUE for masks, FALSE for images. */
    double gamma,         /* Gamma to use in encoding (1 = linear encoding). */    
    double bias,          /* Offset to use in encoding. */                         
    bool_t yup,           /* If TRUE, reverses the indexing of rows. */ 
    bool_t warn,          /* If TRUE, prints "writing {fname}..." to {stderr}. */
    bool_t verbose        /* If TRUE, prints conversion diagnostics to {stderr}. */
  )
  { int32_t i;
    for(i = 0; i < n; i++)
      { float_image_write_pnm_named(fname[i], fim[i], isMask, gamma, bias, yup, warn, verbose); }
  }
