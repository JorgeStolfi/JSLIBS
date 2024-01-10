/* See float_ppm_image.h */
/* Last edited on 2017-06-20 20:47:15 by stolfilocal */ 

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
#include <float_image_from_to_uint16_image.h>
#include <sample_conv.h>
 
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

uint16_image_t *float_image_to_uint16_image
  ( float_image_t *fim,  /* Float image to convert. */
    bool_t isMask,       /* TRUE for masks, FALSE for images. */
    int chns,            /* Number of channels of output image. */
    double lo[],         /* Nominal min float sample for each chosen channel. */
    double hi[],         /* Nominal max float sample for each chosen channel. */
    int ch[],            /* Indices of channels of {fim} to convert. */
    uint16_t maxval, /* Max integer sample value in result image. */
    bool_t yup,          /* If TRUE, reverses the indexing of rows. */
    bool_t verbose       /* If TRUE, prints conversion diagnostics to {stderr}. */
  )
  { /* Get image dimensions: */
    int NX = (int)fim->sz[1];
    int NY = (int)fim->sz[2];
    
    /* Channel counts: */
    int fchns = (int)fim->sz[0]; /* Num channels in float image. */
    int ichns = chns;            /* Num channels in integer image. */
    
    /* Allocate PGM/PPM image: */
    uint16_image_t *pim = uint16_image_new(NX, NY, ichns);
    
    /* Set max sample value in integer image: */
    pim->maxval = maxval;
    
    /* Channel indexing variables: */
    int k; /* Channel of integer image. */
    int c; /* Channel of float image. */
    
    /* Input and output range registers: */
    float vmin[ichns], vmax[ichns];         /* Float pixel range. */
    sample_uint32_t imin[ichns], imax[ichns]; /* Int pixel range. */
    int clo[ichns], chi[ichns];             /* Counts of lo-clipped and hi-clipped pixels. */
    for (k = 0; k < ichns; k++) 
      { clo[k] = chi[k] = 0;
        vmin[k] = +INF;
        vmax[k] = -INF; 
        imin[k] = maxval;
        imax[k] = 0;  
      }
    
    /* Convert pixels, store in {pim}, keep statistics: */
    int x, y;
    for(y = 0; y < NY; y++)
      { int ppmy = (yup ? NY - 1 - y : y);
        uint16_t *prow = pim->smp[ppmy];
        for(x = 0; x < NX; x++)
          { /* Convert float pixel {fpxy[c..c+2]} to integer pixel {ipxy[0..2]}, keep stats: */
            for (k = 0; k < ichns; k++)
              { double lok = (lo == NULL ? 0.0 : lo[k]);
                double hik = (hi == NULL ? 1.0 : hi[k]);
                c = (ch == NULL ? k : ch[k]);
                float v = ((c < 0) || (c >= fchns) ? 0.0f : float_image_get_sample(fim, c, x, y));
                (*prow) = (uint16_t)sample_conv_quantize
                  ( v, maxval, isMask, lok, hik, 
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
              { double lok = (lo == NULL ? 0.0 : lo[k]);
                double hik = (hi == NULL ? 1.0 : hi[k]);
                c = (ch == NULL ? k : ch[k]);
                sample_conv_print_quantize_stats
                  ( c, k, vmin[k], vmax[k], lok, hik, clo[k], chi[k], maxval, imin[k], imax[k]);
              }
          }
      }
    
    return pim;
  }
