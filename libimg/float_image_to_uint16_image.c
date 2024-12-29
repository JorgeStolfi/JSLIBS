/* See {float_image_to_uint16_image.h} */
/* Last edited on 2024-12-26 12:31:44 by stolfi */ 

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
 
#include <float_image_to_uint16_image.h>

uint16_image_t *float_image_to_uint16_image
  ( float_image_t *fim,  /* Float image to convert. */
    bool_t isMask,       /* TRUE for masks, FALSE for images. */
    int32_t chns,            /* Number of channels of output image. */
    double lo[],         /* Nominal min float sample for each chosen channel. */
    double hi[],         /* Nominal max float sample for each chosen channel. */
    int32_t ch[],            /* Indices of channels of {fim} to convert. */
    uint16_t maxval, /* Max integer sample value in result image. */
    bool_t yup,          /* If TRUE, reverses the indexing of rows. */
    bool_t verbose       /* If TRUE, prints conversion diagnostics to {stderr}. */
  )
  { /* Get image dimensions: */
    int32_t NX = (int32_t)fim->sz[1];
    int32_t NY = (int32_t)fim->sz[2];
    
    /* Channel counts: */
    int32_t fchns = (int32_t)fim->sz[0]; /* Num channels in float image. */
    int32_t ichns = chns;            /* Num channels in integer image. */
    
    /* Allocate PGM/PPM image: */
    uint16_image_t *pim = uint16_image_new((uint32_t)NX, (uint32_t)NY, (uint32_t)ichns);
    
    /* Set max sample value in integer image: */
    pim->maxval = maxval;
    
    /* Channel indexing variables: */
    int32_t k; /* Channel of integer image. */
    int32_t c; /* Channel of float image. */
    
    /* Input and output range registers: */
    float vmin[ichns], vmax[ichns];         /* Float pixel range. */
    sample_uint32_t imin[ichns], imax[ichns]; /* Int32_T pixel range. */
    int32_t clo[ichns], chi[ichns];             /* Counts of lo-clipped and hi-clipped pixels. */
    for (k = 0; k < ichns; k++) 
      { clo[k] = chi[k] = 0;
        vmin[k] = +INF;
        vmax[k] = -INF; 
        imin[k] = maxval;
        imax[k] = 0;  
      }
    
    /* Convert pixels, store in {pim}, keep statistics: */
    int32_t x, y;
    for(y = 0; y < NY; y++)
      { int32_t ppmy = (yup ? NY - 1 - y : y);
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
        int32_t NPIX = ((int32_t)NX)*((int32_t)NY);
        fprintf(stderr, "  %d pixels in float image\n", NPIX);
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
