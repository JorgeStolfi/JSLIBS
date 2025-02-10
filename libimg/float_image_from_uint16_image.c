/* See {float_image_from_uint16_image.h} */
/* Last edited on 2025-01-30 08:07:55 by stolfi */ 

#include <limits.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <ix.h>
#include <affirm.h>
#include <bool.h>
 
#include <float_image.h>
#include <sample_conv.h>
 
#include <float_image_from_uint16_image.h>

float_image_t *float_image_from_uint16_image
  ( uint16_image_t *pim, /* PNM image to convert. */
    bool_t isMask,    /* TRUE for masks, FALSE for photos. */
    double lo[],      /* Min float sample value. */
    double hi[],      /* Max float sample value. */
    bool_t yUp,       /* If TRUE, reverses the indexing of rows. */
    bool_t verbose    /* If TRUE, prints conversion diagnostics to {stderr}. */
  )
  { 
    /* Get image dimensions: */
    uint32_t NX = pim->cols;
    uint32_t NY = pim->rows;
    
    /* Channel counts: */
    uint32_t chns = pim->chns; /* Num of channels. */
    
    /* Allocate float image: */
    float_image_t *fim = float_image_new((int32_t)chns, (int32_t)NX, (int32_t)NY);
    
    /* Max sample value in integer image: */
    uint16_t maxval = pim->maxval;
    
    /* Input and output range registers: */
    sample_uint32_t imin[chns], imax[chns]; /* Input range registers. */ 
    float vmin[chns], vmax[chns];         /* Output range registers. */ 
    for (uint32_t c = 0; c < chns; c++) 
      { imin[c] = maxval;
        imax[c] = 0;
        vmin[c] = +INF;
        vmax[c] = -INF;
      }
    
    /* Convert pixels, keep statistics: */
    for(int32_t y = 0; y < NY; y++)
      { int32_t pgmy = (yUp ? (int32_t)NY - 1 - y : y);
        uint16_t *prow = pim->smp[pgmy];
        for(int32_t x = 0; x < NX; x++)
          { for (int32_t c = 0; c < chns; c++)
              { /* Convert int32_t sample {*prow} to float {v}, store, keep stats: */
                uint16_t ismp = (*prow);
                double loc = (lo == NULL ? 0.0 : lo[c]);
                double hic = (hi == NULL ? 1.0 : hi[c]);
                float fsmp = sample_conv_floatize
                  ( ismp, maxval, isMask, loc, hic, &(imin[c]), &(imax[c]), &(vmin[c]), &(vmax[c]) );
                float_image_set_sample(fim, c, x, y, fsmp);
                prow++;
              }
          }
      }
    
    if (verbose) 
      { /* Print statistics: */
        uint32_t NPIX = (NX*NY);
        fprintf(stderr, "  %d pixels in PNM image\n", NPIX);
        if (NPIX > 0)
          { for (int32_t c = 0; c < chns; c++)
              { double loc = (lo == NULL ? 0.0 : lo[c]);
                double hic = (hi == NULL ? 1.0 : hi[c]);
                sample_conv_print_floatize_stats
                  ( c, c, imin[c], imax[c], maxval, loc, hic, vmin[c], vmax[c] );
              }
          }
      }
    
    return fim;
  }

