/* See {float_image_hog.h} */
/* Last edited on 2013-10-21 00:05:39 by stolfilocal */

#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <float_image.h>
#include <affirm.h>

#include <float_image_hog.h>

#define TINY_REL_WEIGHT 1.0e-15
/* If the weight is less than {TINY_REL_WEIGHT} times the total histogram weight, ignore. */

void float_image_hog_collect
  ( float_image_t *DX, 
    int cX,
    float_image_t *DY, 
    int cY,
    float_image_t *M, 
    double noise, 
    bool_t oriented,
    int nh,
    double h[]
  )
  {
    bool_t debug = FALSE;
    
    /* Get the image dimensions: */
    int NX = (int)DX->sz[1]; 
    int NY = (int)DX->sz[2];
    float_image_check_size(DY, -1, NX, NY);
    demand((cX >= 0) && (cX < DX->sz[0]), "invalid {DX} channel");
    demand((cY >= 0) && (cY < DY->sz[0]), "invalid {DY} channel");
    if (M != NULL) { float_image_check_size(M, -1, NX, NY); }
    
    /* Clear histogram: */
    int ih;
    for (ih = 0; ih < nh; ih++) { h[ih] = 0; }
    
    /* Collect histogram: */
    double totw = 0.0; /* Total weight in histogram so far. */
    int x, y;
    for (x = 0; x < NX; x++)
      { for (y = 0; y < NY; y++)
          { /* Compute the horizontal and vertical derivatives in channel {c}: */
            double wM = (M == NULL ? 1.0 : float_image_get_sample(M, 0, x, y));
            if (wM > TINY_REL_WEIGHT*totw)
              { double gx = float_image_get_sample(DX, cX, x, y);
                double gy = float_image_get_sample(DY, cY, x, y);
                /* Compute the gradient modulus squared {gm2}: */
                double gm2 = gx*gx + gy*gy;
                double wg = gm2/(gm2 + noise*noise);
                double w = wg*wM; /* Weight of this sample. */
                if (w > TINY_REL_WEIGHT*totw) 
                  { /* !!! Instead of reducing {wg} for small {gm}, should spread mass {wM} more widely. !!! */
                    /* !!! However that requires finding a good formula for the prob distr of {ga}. !!! */
                    /* Compute the gradient's azimuth: */
                    double ga = atan2(gy,gx);
                    /* Convert {ga} to a fractional position {fh} in {[0_1)} for {[0_PI)} or {[0_2*PI)}: */
                    double fa = ga/(oriented ? 2*M_PI : M_PI);
                    fa = fa - floor(fa);
                    /* Map {fa} to a fractional position {rh} in the domain of {h}: */
                    double rh = nh*fa + 0.5;
                    if (rh >= nh) { rh = rh - nh; }
                    assert((0 <= rh) && (rh < nh));
                    /* Get the bin index and position {sh} in bin: */
                    int ih = (int)floor(rh);
                    assert((0 <= ih) && (ih < nh));
                    double sh = rh - ih;
                    /* Accumulate in adjacent pixels with quadratic interpolation: */
                    double wi = 1 - 2*(sh - 0.5)*(sh - 0.5);
                    h[ih] += w*wi;
                    double wj = 1 - wi;
                    int jh = (sh < 0.5 ? (nh + ih - 1) : ih + 1) % nh;
                    h[jh] += w*wj;
                    if (debug) 
                      { fprintf
                          ( stderr,
                            "G = ( %+9.6f %+9.6f ) a = %+9.4f  w = %8.6f  ih = %3d sh = %5.3f\n", 
                            gx, gy, ga*180/M_PI, w, ih, sh
                            );
                      }
                    totw += w;
                  }
              }
          }
      }
  }
  
  
