/* See float_array_to_uint16_image.h */
/* Last edited on 2017-06-21 22:55:03 by stolfilocal */ 

#define _GNU_SOURCE
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
    int chns,
    int ch[],
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
    int NX = (int)DA->sz[1];  /* Num of columns. */
    int NY = (int)DA->sz[2];  /* Num of rows. */
    
    int iNC = chns;           /* Num channels in integer image. */
    int fNC = (int)DA->sz[0]; /* Num channels in float image. */
    
    /* Allocate PGM/PPM image: */
    uint16_image_t *img = uint16_image_new(NX, NY, iNC);
    
    /* Set max sample value in integer image: */
    demand(maxval + 0u <= uint16_image_MAX_SAMPLE, "invalid maxval");
    img->maxval = (uint16_t)maxval;
    
    /* Input and output range registers: */
    float vmin[iNC], vmax[iNC];         /* Float pixel range. */
    sample_uint32_t imin[iNC], imax[iNC]; /* Int pixel range. */
    int clo[iNC], chi[iNC];             /* Counts of lo-clipped and hi-clipped pixels. */
    for (int ic = 0; ic < iNC; ic++) 
      { clo[ic] = chi[ic] = 0;
        vmin[ic] = +INFINITY;
        vmax[ic] = -INFINITY; 
        imin[ic] = maxval;
        imax[ic] = 0;  
      }
    
    /* Convert pixels, store in {img}, keep statistics: */
    ix_index_t ix[na];
    for(int y = 0; y < NY; y++)
      { /* Fill pixel row {y} of {img}: */
        int fy = (yrev ? NY-1-y : y); /* Row index in float array. */
        uint16_t *prow = img->smp[fy];
        for(int x = 0; x < NX; x++)
          { /* Convert float pixel {fpxy[c..c+2]} to integer pixel {ipxy[0..2]}, keep stats: */
            for (int ic = 0; ic < iNC; ic++)
              { int fc = (ch == NULL ? ic : ch[ic]); /* Channel of float image. */
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
        long int NPIX = ((long int)NX)*((long int)NY);
        fprintf(stderr, "  %ld pixels in float image\n", NPIX);
        if (NPIX > 0)
          { for (int ic = 0; ic < chns; ic++)
              { int fc = (ch == NULL ? ic : ch[ic]);
                double loc = (lo == NULL ? 0.0 : lo[fc]);
                double hic = (hi == NULL ? 1.0 : hi[fc]);
                sample_conv_print_quantize_stats
                  ( fc, ic, vmin[ic], vmax[ic], loc, hic, clo[ic], chi[ic], maxval, imin[ic], imax[ic]);
              }
          }
      }
    
    return img;
  }
