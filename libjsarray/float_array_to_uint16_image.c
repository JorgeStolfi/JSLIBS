/* See float_array_to_uint16_image.h */
/* Last edited on 2024-11-20 07:49:47 by stolfi */ 

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

#include <float_array_to_uint16_image.h>

#define NA  float_array_NAXES
  /* Number of dimensions. */
 
uint16_image_t *float_array_to_uint16_image
  ( float_array_t *A, 
    bool_t isMask,
    int32_t chns,
    int32_t ch[],
    double lo[], 
    double hi[], 
    bool_t yrev, 
    uint16_t maxval, 
    bool_t verbose
  )
  { 
    /* Get indexing descriptor {DA}: */
    ix_descr_t *DA = &(A->ds);
    
    /* Get integer image dimensions: */
    ix_dim_t na = DA->na;
    assert(na == 3);
    int32_t NX = (int32_t)DA->sz[1];  /* Num of columns. */
    int32_t NY = (int32_t)DA->sz[2];  /* Num of rows. */
    
    int32_t iNC = chns;           /* Num channels in integer image. */
    int32_t fNC = (int32_t)DA->sz[0]; /* Num channels in float image. */
    
    /* Allocate PGM/PPM image: */
    uint16_image_t *img = uint16_image_new(NX, NY, iNC);
    
    /* Set max sample value in integer image: */
    demand(maxval + 0u <= uint16_image_MAX_SAMPLE, "invalid maxval");
    img->maxval = (uint16_t)maxval;
    
    /* Input and output range registers: */
    float vmin[iNC], vmax[iNC];         /* Float pixel range. */
    sample_uint32_t imin[iNC], imax[iNC]; /* Int32_T pixel range. */
    int32_t clo[iNC], chi[iNC];             /* Counts of lo-clipped and hi-clipped pixels. */
    for (uint32_t ic = 0;  ic < iNC; ic++) 
      { clo[ic] = chi[ic] = 0;
        vmin[ic] = +INFINITY;
        vmax[ic] = -INFINITY; 
        imin[ic] = maxval;
        imax[ic] = 0;  
      }
    
    /* Convert pixels, store in {img}, keep statistics: */
    ix_index_t ix[na];
    for (uint32_t y = 0;  y < NY; y++)
      { /* Fill pixel row {y} of {img}: */
        int32_t fy = (yrev ? NY-1-y : y); /* Row index in float array. */
        uint16_t *prow = img->smp[fy];
        for (uint32_t x = 0;  x < NX; x++)
          { /* Convert float pixel {fpxy[c..c+2]} to integer pixel {ipxy[0..2]}, keep stats: */
            for (uint32_t ic = 0;  ic < iNC; ic++)
              { int32_t fc = (ch == NULL ? ic : ch[ic]); /* Channel of float image. */
                ix[0] = fc; ix[1] = x; ix[2] = y; ix[3] = ix[4] = ix[5] = 0;
                float v = ((fc < 0) || (fc >= fNC) ? 0.0f : float_array_get_elem(A, ix));
                double loc = (lo == NULL ? 0.0 : lo[fc]);
                double hic = (hi == NULL ? 1.0 : hi[fc]);
                (*prow) = (uint16_t)sample_conv_quantize
                  ( v, maxval, isMask, loc, hic, 
                    &(vmin[ic]), &(vmax[ic]), &(clo[ic]), &(chi[ic]), &(imin[ic]), &(imax[ic])
                  );
                prow++;
              }
          }
      }
    
    if (verbose)
      { /* Print statistics: */
        int64_t NPIX = ((int64_t)NX)*((int64_t)NY);
        fprintf(stderr, "  %ld pixels in float image\n", NPIX);
        if (NPIX > 0)
          { for (uint32_t ic = 0;  ic < chns; ic++)
              { int32_t fc = (ch == NULL ? ic : ch[ic]);
                double loc = (lo == NULL ? 0.0 : lo[fc]);
                double hic = (hi == NULL ? 1.0 : hi[fc]);
                sample_conv_print_quantize_stats
                  ( fc, ic, vmin[ic], vmax[ic], loc, hic, clo[ic], chi[ic], maxval, imin[ic], imax[ic]);
              }
          }
      }
    
    return img;
  }
