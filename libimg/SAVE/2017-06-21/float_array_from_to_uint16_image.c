/* See float_array_from_to_uint16_image.h */
/* Last edited on 2017-06-21 01:22:22 by stolfilocal */ 

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
#include <jspnm.h>
#include <uint16_image.h>
 
#include <float_array.h>
#include <float_array_from_to_uint16_image.h>
#include <sample_conv.h>

#define NA  float_array_NAXES
  /* Number of dimensions. */
 
float_array_t *float_array_from_uint16_image
  ( uint16_image_t *img, 
    bool_t isMask,
    double lo[], 
    double hi[], 
    bool_t verbose
  )
  { 
    /* Get image dimensions: */
    int NX = img->cols;
    int NY = img->rows;
    
    /* Channel counts: */
    int NC = img->chns; /* Num of channels. */
    
    /* Allocate float image: */
    ix_dim_t na = 3;
    ix_size_t sz[na];
    sz[0] = NC; sz[1] = NX; sz[2] = NY; 
    float_array_t AN = float_array_new(na, sz, ix_order_L);
    float_array_t *A = &AN;
    
    /* Max sample value in integer image: */
    uint16_t maxval = img->maxval;
    
    /* Input and output range registers: */
    sample_uint32_t imin[NC], imax[NC]; /* Input range registers. */ 
    float vmin[NC], vmax[NC];    /* Output range registers. */ 
    int c; /* Channel index. */
    for (c = 0; c < NC; c++) 
      { imin[c] = maxval;
        imax[c] = 0;
        vmin[c] = +INFINITY;
        vmax[c] = -INFINITY;
      }
    
    /* Convert pixels, keep statistics: */
    ix_index_t ix[na];
    int x, y;
    for(y = 0; y < NY; y++)
      { int pgmy = NY - 1 - y;
        uint16_t *prow = img->smp[pgmy];
        for(x = 0; x < NX; x++)
          { for (c = 0; c < NC; c++)
              { /* Convert int sample {*prow} to float {v}, store, keep stats: */
                sample_uint32_t ismp = (*prow);
                double loc = (lo == NULL ? 0.0 : lo[c]);
                double hic = (hi == NULL ? 1.0 : hi[c]);
                float fsmp = sample_conv_floatize
                  ( ismp, maxval, isMask, loc, hic, &(imin[c]), &(imax[c]), &(vmin[c]), &(vmax[c]) );
                ix[0] = c; ix[1] = x; ix[2] = y; ix[3] = ix[4] = ix[5] = 0;
                float_array_set_elem(A, ix, fsmp);
                prow++;
              }
          }
      }
    
    if (verbose) 
      { /* Print statistics: */
        long int NPIX = ((long int)NX)*((long int)NY);
        fprintf(stderr, "  %ld pixels in PNM image\n", NPIX);
        if (NPIX > 0)
          { for (c = 0; c < NC; c++)
              { double loc = (lo == NULL ? 0.0 : lo[c]);
                double hic = (hi == NULL ? 1.0 : hi[c]);
                sample_conv_print_floatize_stats
                  ( c, c, imin[c], imax[c], maxval, loc, hic, vmin[c], vmax[c] );
              }
          }
      }
    
    return A;
  }

uint16_image_t *float_array_to_uint16_image
  ( float_array_t *A, 
    bool_t isMask,
    int chns,
    int ch[],
    double lo[], 
    double hi[], 
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
    uint16_image_t *iim = uint16_image_new(NX, NY, iNC);
    
    /* Set max sample value in integer image: */
    demand(maxval + 0u <= uint16_image_MAX_SAMPLE, "invalid maxval");
    iim->maxval = (uint16_t)maxval;
    
    /* Channel indexing variables: */
    int k; /* Channel of integer image. */
    int c; /* Channel of float image. */
    
    /* Input and output range registers: */
    float vmin[iNC], vmax[iNC];         /* Float pixel range. */
    sample_uint32_t imin[iNC], imax[iNC]; /* Int pixel range. */
    int clo[iNC], chi[iNC];             /* Counts of lo-clipped and hi-clipped pixels. */
    for (k = 0; k < iNC; k++) 
      { clo[k] = chi[k] = 0;
        vmin[k] = +INFINITY;
        vmax[k] = -INFINITY; 
        imin[k] = maxval;
        imax[k] = 0;  
      }
    
    /* Convert pixels, store in {iim}, keep statistics: */
    ix_index_t ix[na];
    int x, y;
    for(y = 0; y < NY; y++)
      { int ppmy = NY-1-y;
        uint16_t *prow = iim->smp[ppmy];
        for(x = 0; x < NX; x++)
          { /* Convert float pixel {fpxy[c..c+2]} to integer pixel {ipxy[0..2]}, keep stats: */
            for (k = 0; k < iNC; k++)
              { c = (ch == NULL ? k : ch[k]);
                ix[0] = c; ix[1] = x; ix[2] = y; ix[3] = ix[4] = ix[5] = 0;
                float v = ((c < 0) || (c >= fNC) ? 0.0f : float_array_get_elem(A, ix));
                double loc = (lo == NULL ? 0.0 : lo[c]);
                double hic = (hi == NULL ? 1.0 : hi[c]);
                (*prow) = (uint16_t)sample_conv_quantize
                  ( v, maxval, isMask, loc, hic, 
                    &(vmin[k]), &(vmax[k]), &(clo[k]), &(chi[k]), &(imin[k]), &(imax[k])
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
          { for (k = 0; k < chns; k++)
              { c = (ch == NULL ? k : ch[k]);
                double loc = (lo == NULL ? 0.0 : lo[c]);
                double hic = (hi == NULL ? 1.0 : hi[c]);
                sample_conv_print_quantize_stats
                  ( c, k, vmin[k], vmax[k], loc, hic, clo[k], chi[k], maxval, imin[k], imax[k]);
              }
          }
      }
    
    return iim;
  }
