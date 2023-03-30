/* See float_array_from_uint16_image.h */
/* Last edited on 2023-03-19 15:25:53 by stolfi */ 

#define _GNU_SOURCE
#include <stdint.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <ix.h>
#include <affirm.h>
#include <bool.h>
#include <jspnm.h>
#include <uint16_image.h>
#include <float_array.h>
#include <sample_conv.h>

#include <float_array_from_uint16_image.h>

#define NA  float_array_NAXES
  /* Number of dimensions. */
 
float_array_t *float_array_from_uint16_image
  ( uint16_image_t *img, 
    bool_t isMask,
    double lo[], 
    double hi[], 
    bool_t yrev,
    bool_t verbose
  )
  { 
    /* Get image dimensions: */
    int32_t NX = img->cols;
    int32_t NY = img->rows;
    
    /* Channel counts: */
    int32_t NC = img->chns; /* Num of channels. */
    
    /* Allocate float image: */
    ix_dim_t na = 3;
    ix_size_t sz[na];
    sz[0] = NC; sz[1] = NX; sz[2] = NY; 
    float_array_t AN = float_array_new(na, sz);
    float_array_t *A = &AN;
    
    /* Max sample value in integer image: */
    uint16_t maxval = img->maxval;
    
    /* Input and output range registers: */
    sample_uint32_t imin[NC], imax[NC]; /* Input range registers. */ 
    float vmin[NC], vmax[NC];    /* Output range registers. */ 
    for (int32_t c = 0; c < NC; c++) 
      { imin[c] = maxval;
        imax[c] = 0;
        vmin[c] = +INFINITY;
        vmax[c] = -INFINITY;
      }
    
    /* Convert pixels, keep statistics: */
    ix_index_t ix[na];
    for(int32_t fy = 0; fy < NY; fy++)
      { /* Fill array plane {ix[2]==fy}: */
        int32_t iy = (yrev ? NY - 1 - fy : fy); /* Row index in image array. */
        uint16_t *prow = img->smp[iy];
        for(int32_t x = 0; x < NX; x++)
          { for (int32_t c = 0; c < NC; c++)
              { /* Convert int32_t sample {*prow} to float {v}, store, keep stats: */
                sample_uint32_t ismp = (*prow);
                double loc = (lo == NULL ? 0.0 : lo[c]);
                double hic = (hi == NULL ? 1.0 : hi[c]);
                float fsmp = sample_conv_floatize
                  ( ismp, maxval, isMask, loc, hic, &(imin[c]), &(imax[c]), &(vmin[c]), &(vmax[c]) );
                ix[0] = c; ix[1] = x; ix[2] = fy; ix[3] = ix[4] = ix[5] = 0;
                float_array_set_elem(A, ix, fsmp);
                prow++;
              }
          }
      }
    
    if (verbose) 
      { /* Print statistics: */
        int64_t NPIX = ((int64_t)NX)*((int64_t)NY);
        fprintf(stderr, "  %ld pixels in PNM image\n", NPIX);
        if (NPIX > 0)
          { for (int32_t c = 0; c < NC; c++)
              { double loc = (lo == NULL ? 0.0 : lo[c]);
                double hic = (hi == NULL ? 1.0 : hi[c]);
                sample_conv_print_floatize_stats
                  ( c, c, imin[c], imax[c], maxval, loc, hic, vmin[c], vmax[c] );
              }
          }
      }
    
    return A;
  }
