/* See {float_image_io_pnm.h}. */
/* Last edited on 2017-06-21 05:24:51 by stolfilocal */

#include "uint16_image.h"
#include "float_image_from_to_uint16_image.h"
#include <assert.h>
#include <affirm.h>
#include <bool.h>

#include "float_image_io_pnm.h"

float_image_t *float_image_from_to_uint16_image_read
  ( char *fname,    /* PPM/PGM/PBM file name (with extension). */
    bool_t isMask,  /* TRUE for masks, FALSE for images. */
    double gamma,   /* Gamma to use in decoding (0 = linear decoding). */
    double bias,    /* Offset to use in decoding. */
    bool_t yup,     /* If TRUE, reverses the indexing of rows. */
    bool_t warn,    /* If TRUE, prints "reading {fname}..." to {stderr}. */
    bool_t verbose  /* If TRUE, prints conversion diagnostics to {stderr}. */
  )
  { uint16_image_t *pim = uint16_image_read(fname, warn);
    float_image_t *fim = float_image_from_uint16_image(pim, isMask, NULL, NULL, yup, verbose);
    if (gamma != 1)
      { int c;
        for (c = 0; c < fim->sz[0]; c++)
          { float_image_apply_gamma(fim, c, gamma, bias); }
      }
    uint16_image_free(pim);
    return fim;
  }

void float_image_from_to_uint16_image_write
  ( char *fname,        /* PPM/PGM/PBM file name (with extension). */            
    float_image_t *fim, /* Image to write. */
    bool_t isMask,      /* TRUE for masks, FALSE for images. */
    double gamma,       /* Gamma to use in encoding (1 = linear encoding). */    
    double bias,        /* Offset to use in encoding. */                         
    bool_t yup,         /* If TRUE, reverses the indexing of rows. */            
    bool_t warn,        /* If TRUE, prints "writing {fname}..." to {stderr}. */
    bool_t verbose      /* If TRUE, prints conversion diagnostics to {stderr}. */
  )
  { int nc = (int)fim->sz[0];
    assert((nc == 1) || (nc == 3));
    float_image_t * gim = float_image_copy(fim);
    int c;
    for (c = 0; c < fim->sz[0]; c++)
      { float_image_apply_gamma(gim, c, 1.0/gamma, bias); }
    uint16_t maxval = PNM_FILE_MAX_MAXVAL;
    uint16_image_t *pim = 
      float_image_to_uint16_image(gim, isMask, nc, NULL, NULL, NULL, maxval, yup, verbose);
    float_image_free(gim);
    uint16_image_write(fname, pim, FALSE, warn);
    uint16_image_free(pim);
  }

float_image_t **float_image_from_to_uint16_image_list_read
  ( int n,          /* Number of images to read. */
    char *fname[],  /* PPM/PGM/PBM file names (with extensions). */
    bool_t isMask,  /* TRUE for masks, FALSE for images. */
    double gamma,   /* Gamma to use in decoding (1 = linear decoding). */
    double bias,    /* Offset to use in decoding. */
    bool_t yup,     /* If TRUE, reverses the indexing of rows. */
    bool_t warn,    /* If TRUE, prints "reading {fname}..." to {stderr}. */
    bool_t verbose  /* If TRUE, prints conversion diagnostics to {stderr}. */
  )
  { float_image_t **fim = notnull(malloc(n * sizeof(float_image_t *)), "no mem");
    int i;
    for(i = 0; i < n; i++)
      { fim[i] = float_image_from_to_uint16_image_read(fname[i], isMask, gamma, bias, yup, warn, verbose); }
    return fim;
  }

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
  )
  { int i;
    for(i = 0; i < n; i++)
      { float_image_from_to_uint16_image_write(fname[i], fim[i], isMask, gamma, bias, yup, warn, verbose); }
  }
