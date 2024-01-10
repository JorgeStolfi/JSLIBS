/* See float_array_pnm.h */
/* Last edited on 2009-08-31 11:29:14 by stolfi */ 

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
#include <jspnm_image.h>
 
#include <float_array.h>
#include <float_array_pnm.h>
#include <sample_conv.h>

#define NA  float_array_NAXES
  /* Number of dimensions. */
 
float_array_t float_array_from_pnm_image
  ( pnm_image_t *img, 
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
    float_array_t A = float_array_new(na, sz, ix_order_L);
    
    /* Max sample value in integer image: */
    pnm_sample_t maxval = img->maxval;
    
    /* Input and output range registers: */
    sample_uint_t imin[NC], imax[NC]; /* Input range registers. */ 
    float vmin[NC], vmax[NC];         /* Output range registers. */ 
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
        pnm_sample_t *prow = img->smp[pgmy];
        for(x = 0; x < NX; x++)
          { for (c = 0; c < NC; c++)
              { /* Convert int sample {*prow} to float {v}, store, keep stats: */
                pnm_sample_t ismp = (*prow);
                double loc = (lo == NULL ? 0.0 : lo[c]);
                double hic = (hi == NULL ? 1.0 : hi[c]);
                float fsmp = sample_conv_floatize
                  ( ismp, maxval, loc, hic, &(imin[c]), &(imax[c]), &(vmin[c]), &(vmax[c]) );
                ix[0] = c; ix[1] = x; ix[2] = y; ix[3] = ix[4] = ix[5] = 0;
                float_array_set_element(&A, ix, fsmp);
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

pnm_image_t *float_array_to_pnm_image
  ( float_array_t *A, 
    int chns,
    int ch[],
    double lo[], 
    double hi[], 
    pnm_sample_t maxval, 
    bool_t verbose
  )
  { 
    /* Get indexing descriptor {D}: */
    ix_descr_t *D = &(A->ds);
    
    /* Get integer image dimensions: */
    ix_dim_t na = D->na;
    assert(na == 3);
    int NX = D->sz[1];  /* Num of columns. */
    int NY = D->sz[2];  /* Num of rows. */
    
    int iNC = chns;     /* Num channels in integer image. */
    int fNC = D->sz[0]; /* Num channels in float image. */
    
    /* Allocate PGM/PPM image: */
    pnm_image_t *iim = pnm_image_new(NX, NY, iNC);
    
    /* Set max sample value in integer image: */
    iim->maxval = maxval;
    
    /* Channel indexing variables: */
    int k; /* Channel of integer image. */
    int c; /* Channel of float image. */
    
    /* Input and output range registers: */
    float vmin[iNC], vmax[iNC];         /* Float pixel range. */
    sample_uint_t imin[iNC], imax[iNC]; /* Int pixel range. */
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
        pnm_sample_t *prow = iim->smp[ppmy];
        for(x = 0; x < NX; x++)
          { /* Convert float pixel {fpxy[c..c+2]} to integer pixel {ipxy[0..2]}, keep stats: */
            for (k = 0; k < iNC; k++)
              { c = (ch == NULL ? k : ch[k]);
                ix[0] = c; ix[1] = x; ix[2] = y; ix[3] = ix[4] = ix[5] = 0;
                float v = ((c < 0) || (c >= fNC) ? 0.0 : float_array_get_element(A, ix));
                double loc = (lo == NULL ? 0.0 : lo[c]);
                double hic = (hi == NULL ? 1.0 : hi[c]);
                (*prow) = sample_conv_quantize
                  ( v, maxval, loc, hic, 
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
