/* See {float_image_read_pnm.h}. */
/* Last edited on 2025-01-30 08:07:35 by stolfi */

#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <float_image.h>
#include <uint16_image.h>
#include <uint16_image_read_pnm.h>
#include <float_image_from_uint16_image.h>

#include <float_image_read_pnm.h>

float_image_t *float_image_read_pnm_named
  ( char *fname,    /* PPM/PGM/PBM file name (with extension). */
    bool_t isMask,  /* TRUE for masks, FALSE for images. */
    double gamma,   /* Gamma to use in decoding (0 = linear decoding). */
    double bias,    /* Offset to use in decoding. */
    bool_t yUp,     /* If TRUE, reverses the indexing of rows. */
    bool_t warn,    /* If TRUE, prints "reading {fname}..." to {stderr}. */
    bool_t verbose  /* If TRUE, prints conversion diagnostics to {stderr}. */
  )
  { uint16_image_t *pim = uint16_image_read_pnm_named(fname, warn);
    float_image_t *fim = float_image_from_uint16_image(pim, isMask, NULL, NULL, yUp, verbose);
    if ((! isnan(gamma)) && (gamma >= 0) && (gamma != 1))
      { int32_t c;
        for (c = 0; c < fim->sz[0]; c++)
          { float_image_apply_gamma(fim, c, gamma, bias); }
      }
    uint16_image_free(pim);
    return fim;
  }

float_image_t **float_image_read_pnm_named_list
  ( int32_t n,          /* Number of images to read. */
    char *fname[],  /* PPM/PGM/PBM file names (with extensions). */
    bool_t isMask,  /* TRUE for masks, FALSE for images. */
    double gamma,   /* Gamma to use in decoding (1 = linear decoding). */
    double bias,    /* Offset to use in decoding. */
    bool_t yUp,     /* If TRUE, reverses the indexing of rows. */
    bool_t warn,    /* If TRUE, prints "reading {fname}..." to {stderr}. */
    bool_t verbose  /* If TRUE, prints conversion diagnostics to {stderr}. */
  )
  { float_image_t **fim = talloc(n, float_image_t *); 
    int32_t i;
    for(i = 0; i < n; i++)
      { fim[i] = float_image_read_pnm_named(fname[i], isMask, gamma, bias, yUp, warn, verbose); }
    return fim;
  }

