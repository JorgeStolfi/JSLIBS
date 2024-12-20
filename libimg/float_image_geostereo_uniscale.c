/* See {float_image_geostereo.h}. */
/* Last edited on 2024-12-05 00:48:30 by stolfi */

#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
 
#include <bool.h>
#include <affirm.h>
#include <wt_table.h>
#include <wt_table_binomial.h>
#include <float_image.h>
#include <float_image_geostereo.h>

#include <float_image_geostereo_uniscale.h>

#define HDEBUG 300
#define VDEBUG 300
  /* Print debugging information for this pixel. */

void float_image_geostereo_uniscale
  ( float_image_t *f1,  /* Image 1. */
    float_image_t *f2,  /* Image 2. */
    int32_t nwx,         /* Window width. */
    int32_t nwy,         /* Window height. */
    double dmin,         /* Minimum signed displacement (pixels). */
    double dmax,         /* Maximum signed displacement (pixels). */
    int32_t ncands,      /* Number of candidates to keep. */
    float_image_t **fd,  /* (OUT) Dispmap image. */
    float_image_t **fs   /* (OUT) Scoremap image. */
  )
  {
    /* Get and check image and window sizes: */
    int32_t NC, NX, NY;
    float_image_get_size(f1 ,&NC, &NX, &NY);
    float_image_check_size(f2, NC, NX, NY);
    demand((nwx % 2) == 1, "window width must be odd");
    demand((nwy % 2) == 1, "window height must be odd");
    int32_t npix = nwx*nwy; /* Number of pixels in window. */
    int32_t nsmp = NC*npix; /* Number of samples in window. */
    
    /* Queue of best candidates for each pixel: */
    double dbest[ncands];
    double sbest[ncands];
    
    /* Window pixel weight table: */
    double *wt = talloc(npix, double);
    double wtx[nwx]; wt_table_binomial_fill((uint32_t)nwx, wtx, NULL);
    wt_table_normalize_sum((uint32_t)nwx, wtx);
    double wty[nwy]; wt_table_binomial_fill((uint32_t)nwy, wty, NULL);
    wt_table_normalize_sum((uint32_t)nwy, wty);
    for (int32_t iy = 0;  iy < nwy; iy++)
      { for (int32_t ix = 0;  ix < nwx; ix++)
          { wt[ix + nwx*iy] = wtx[ix]*wty[iy]; }
      }
    
    /* Work areas for window samples: */
    double *smp1 = talloc(nsmp, double);
    double *smp2 = talloc(nsmp, double);

    /* Allocate dispmap and scoremap: */
    (*fd) = float_image_new(ncands, NX, NY);
    (*fs) = float_image_new(ncands, NX, NY);

    /* Compute dispmap/scoremap: */
    for (int32_t y = 0;  y < NY; y++)
      { for (int32_t x = 0;  x < NX; x++)
         { bool_t debug = ((x == HDEBUG) && (y == VDEBUG));
           float_image_geostereo_single_pixel_best
             ( f1, f2, 
               x, y, 
               nwx, nwy, wt, 
               dmin, dmax, 
               ncands, dbest, sbest, 
               debug, 
               smp1, smp2
             );

           /* Store best candidates in heigh map and score map: */
           for (int32_t rk = 0;  rk < ncands; rk++)
             { float_image_set_sample((*fd), rk, x, y, (float)dbest[rk]);
               float_image_set_sample((*fs), rk, x, y, (float)sbest[rk]);
             }
         }
      }
      
    /* Free auxiliary storage: */
    free(wt);
    free(smp1);
    free(smp2);
  }
