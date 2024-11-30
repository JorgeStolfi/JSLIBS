/* See {float_image_noise.h}. */

/* Last edited on 2024-10-26 04:27:24 by stolfi */ 
/* Created on 2009-06-02 by J. Stolfi, UNICAMP */

#define float_image_test_C_COPYRIGHT \
  "Copyright © 2024  by the State University of Campinas (UNICAMP)"

#include <math.h>
#include <stdint.h>
#include <assert.h>

#include <affirm.h>
#include <jsrandom.h>
#include <gauss_table.h>

#include <float_image.h>
#include <float_image_hartley.h>

#include <float_image_noise.h>

void float_image_noise_fill_filter_table(double ff, int32_t NW, double wt[]);
  /* Fills the filter table {wt[0..NW-1]} with a folded Gaussian-like 1D
    weight function with deviation {ff}. The {ff} parameter must be
    finite and non-negative. If {ff} is 0, the table will be all 1.0. */

float_image_t *float_image_noise
  ( int32_t NC,
    int32_t NX,
    int32_t NY,
    double fxFilter,
    double fyFilter,
    bool_t complement
  )
  {
    /* Allocate the transform image: */
    float_image_t *img = float_image_new(NC, NX,NY);
    
    /* Set all componenents to unit aplitude and random phase: */
    for (uint32_t ic = 0;  ic < NC; ic++)
      { for (int32_t ix = 0; ix < NX; ix ++)
          { int32_t jx = (NX - ix) % NX;
            for (int32_t iy = 0; iy < NY; iy++ )
              { int32_t jy = (NY - iy) % NY;
                if ((ix == 0) && (iy == 0))
                  { /* Constant term. Omit it: */
                    float_image_set_sample(img, ic, ix, iy, 0.0);
                  }
                else if ((ix == jx) && (iy == jy))
                  { /* A self-conjugated frequency. Set it to unit amplitude: */
                    float_image_set_sample(img, ic, ix, iy, 1.0);
                  }
                else if ((ix < jx) || ((ix == jx) && (iy <= jy)))
                  { assert ((ix != jx) || (iy != jy));
                    double ang = drandom()*2*M_PI;
                    double sa = sin(ang), ca = cos(ang);
                    float_image_set_sample(img, ic, ix, iy, (float)ca);
                    float_image_set_sample(img, ic, jx, jy, (float)sa);
                  }
              }
          }
      }
      
    /* Apply filters: */
    if (fxFilter == +INF) { fxFilter = 0.0; }
    if (fyFilter == +INF) { fyFilter = 0.0; }
    if ((fxFilter != 0) || (fyFilter != 0) || complement)
      { double wtx[NX]; float_image_noise_fill_filter_table(fxFilter, NX, wtx);
        double wty[NY]; float_image_noise_fill_filter_table(fyFilter, NY, wty);
        for (int32_t ix = 0; ix < NX; ix ++)
          { for (int32_t iy = 0; iy < NY; iy++ )
              { double wxy = wtx[ix]*wty[iy];
                if (complement) { wxy = 1 - wxy; }
                for (uint32_t ic = 0;  ic < NC; ic++)
                  { double smp = float_image_get_sample(img, ic, ix, iy);
                    float_image_set_sample(img, ic, ix, iy, (float)(wxy*smp));
                  }
              }
          }
      }
    
    /* Inverse transform: */
    float_image_hartley_transform(img, img);
    return img;
  }
    
void float_image_noise_fill_filter_table(double ff, int32_t NW, double wt[])
  { 
    bool_t debug = TRUE;
    
    demand(isfinite(ff) && (ff >= 0), "invalid filter freq parameter");
    if (ff == 0)
      { for (uint32_t iw = 0;  iw < NW; iw++) { wt[iw] = 1.0; } }
    else
      { double mag = -INF;
        for (uint32_t iw = 0;  iw < NW; iw++)
          { wt[iw] = gauss_table_folded_bell((double)iw, ff, NW);
            assert(wt[iw] >= 0);
            mag = fmax(mag, wt[iw]);
          }
        for (uint32_t iw = 0;  iw < NW; iw++) 
          { wt[iw] /= mag;
            if (debug) { fprintf(stderr, "  wt[%4d] = %10.8f\n", iw, wt[iw]); }
          }
      }
  }
