/* See {float_image_from_uint16_image.h} */
/* Last edited on 2017-06-22 17:43:19 by stolfilocal */ 

#define _GNU_SOURCE
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <indexing.h>
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
    bool_t yup,       /* If TRUE, reverses the indexing of rows. */
    bool_t verbose    /* If TRUE, prints conversion diagnostics to {stderr}. */
  )
  { 
    /* Get image dimensions: */
    int NX = pim->cols;
    int NY = pim->rows;
    
    /* Channel counts: */
    int chns = pim->chns; /* Num of channels. */
    
    /* Allocate float image: */
    float_image_t *fim = float_image_new(chns, NX, NY);
    
    /* Max sample value in integer image: */
    uint16_t maxval = pim->maxval;
    
    /* Input and output range registers: */
    sample_uint32_t imin[chns], imax[chns]; /* Input range registers. */ 
    float vmin[chns], vmax[chns];         /* Output range registers. */ 
    int c; /* Channel index. */
    for (c = 0; c < chns; c++) 
      { imin[c] = maxval;
        imax[c] = 0;
        vmin[c] = +INF;
        vmax[c] = -INF;
      }
    
    /* Convert pixels, keep statistics: */
    int x, y;
    for(y = 0; y < NY; y++)
      { int pgmy = (yup ? NY - 1 - y : y);
        uint16_t *prow = pim->smp[pgmy];
        for(x = 0; x < NX; x++)
          { for (c = 0; c < chns; c++)
              { /* Convert int sample {*prow} to float {v}, store, keep stats: */
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
        long int NPIX = ((long int)NX)*((long int)NY);
        fprintf(stderr, "  %ld pixels in PNM image\n", NPIX);
        if (NPIX > 0)
          { for (c = 0; c < chns; c++)
              { double loc = (lo == NULL ? 0.0 : lo[c]);
                double hic = (hi == NULL ? 1.0 : hi[c]);
                sample_conv_print_floatize_stats
                  ( c, c, imin[c], imax[c], maxval, loc, hic, vmin[c], vmax[c] );
              }
          }
      }
    
    return fim;
  }

