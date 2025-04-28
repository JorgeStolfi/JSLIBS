/* See {multifok_weight.h}. */
/* Last edited on 2025-04-13 15:24:48 by stolfi */

#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <r3.h>
#include <float_image.h>

#include <multifok_weight.h>

float_image_t *multifok_weight_from_height_dev(float_image_t *hDev, double dhMax)
  {
    demand(hDev != NULL, "{hDev} must be non-null");
    
    int32_t NC, NX, NY;
    float_image_get_size(hDev, &NC, &NX, &NY);
      
    demand(NC >= 1, "invalid {hDev} channel count");
        
    auto double compute_weight(int32_t x, int32_t y);
      /* Return the weight for pixel {[x,y]} (either 0 or 1) 
        considering {hDev[0,x,y]} and {hdev[1,x,y]} if 
        it exists. */

    float_image_t *hWht = float_image_new(1, NX, NY);
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { double w = compute_weight(x, y);
            float_image_set_sample(hWht, 0, x, y, (float)w);
          }
      }
      
    return hWht;
    
    double compute_weight(int32_t x, int32_t y)
      { /* Check height deviation: */
        float dh = float_image_get_sample(hDev, 0, x, y);
        if (! isfinite(dh)) { return 0.0; }
        demand(dh >= 0, "deviation cannot be negative");
        if (dh > dhMax) { return 0.0; }
        
        if (NC >= 2)
          { /* Get the weight from {hDev} too: */
            float wh = float_image_get_sample(hDev, 1, x, y);
            if (! isnan(wh)) 
              { demand(isfinite(wh) && (wh >= 0), "invalid weight value in height deviation map");
                return wh;
              }
          }
        return 1.0;
      }
  }

float_image_t *multifok_weight_from_normal_grad(float_image_t *sNrm, double dnMax)
  {
    demand(sNrm != NULL, "{sNrm} must not be null");
    
    int32_t NC, NX, NY;
    float_image_get_size(sNrm, &NC, &NX, &NY);
      
    demand(NC >= 3, "invalid {sNrm} channel count");
        
    auto double compute_weight(int32_t x, int32_t y);
      /* Return the weight for pixel {[x,y]} (either 0 or 1) considering {sNrm} alone.
        In particular, if {sNrm} is null, returns 1. */
        
    auto bool_t get_normal(int32_t x, int32_t y, r3_t *u_P);
      /* Tries to gets the normal vector {u=sNrm[0..2,x,y]} and normalize
        it to unit length. If the pixel {[x,y]} does not exist or {u} is
        any coordinate of {u} is {NAN} or the norm of {u} is too close
        to zero, returns {FALSE} and {*u_P} is undefined. Othwerwise
        stores the normalized {u} into {*u_P} and returns {TRUE} */

    float_image_t *nWht = float_image_new(1, NX, NY);
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { double w = compute_weight(x, y);
            float_image_set_sample(nWht, 0, x, y, (float)w);
          }
      }
      
    return nWht;
    
    double compute_weight(int32_t x, int32_t y)
      { /* Get normal at pixel: */
        r3_t u0;  
        bool_t u0_ok = get_normal(x, y, &u0);
        if (! u0_ok) { return 0.0; }

        /* Compare with adjacent pixels: */
        for (int32_t ey = -1; ey <= +1; ey++)
          { for (int32_t ex = -1; ex <= +1; ex++)
              { int32_t x1 = x + ex, y1 = y + ey;
                r3_t u1;
                bool_t u1_ok = get_normal(x1, y1, &u1);
                if (u1_ok)
                  { double dn = r3_dist(&u0, &u1);
                    if (dn > dnMax) { return 0.0; }
                  }
              }
          }
        if (NC >= 4)
          { /* Get the weight from {sNrm} too: */
            float wn = float_image_get_sample(sNrm, 1, x, y);
            if (! isnan(wn))
              { demand(isfinite(wn) && (wn >= 0), "invalid weight value in normal map");
                return wn;
              }
          }
        return 1.0;
      }
        
    bool_t get_normal( int32_t x, int32_t y, r3_t *u_P)
      { if ((x < 0) || (x >= NX) || (y < 0) || (y >= NY)) { return FALSE; }
        r3_t u;
        for (int32_t c = 0; c < 3; c++) 
          { u.c[c] = float_image_get_sample(sNrm, c, x, y);
            if (! isfinite(u.c[c])) { return FALSE; } 
            /* No need to check for {±INF}, will get to it eventually. */
          }
        /* Try to normalize, just in case: */
        double u_mag = r3_dir(&u, &u);
        if (u_mag < 1.0e-4) { return FALSE; }
        (*u_P) = u;
        return TRUE;
      }
  }

void multifok_weight_multiply(float_image_t *dst, int32_t cd, float_image_t *src, int32_t cs)
  {
    bool_t debug = TRUE;
    
    int32_t NC_d, NC_s, NX, NY;
    float_image_get_size(dst, &NC_d, &NX, &NY);
    float_image_check_size(src, -1, NX, NY, "images {src,dst} have different sizes");
    if ((cd < 0) || (cd >= NC_d)) { return; }
    
    NC_s = (int32_t)(src->sz[0]);
    if ((cs < 0) || (cs >= NC_s)) { return; }
    
    if (debug) { fprintf(stderr, "entered %s ...\n", __FUNCTION__); }

    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { float *wdp = float_image_get_sample_address(dst, cd, x, y);
            double wd = (*wdp);
            if (! isnan(wd))
              { demand(isfinite(wd) && (wd >= 0), "invalid weight value in {dst} map"); };
            double ws = float_image_get_sample(src, cs, x, y);
            if (! isnan(ws))
              { demand(isfinite(ws) && (ws >= 0), "invalid weight value in {src} map"); };
            double w;
            if (isnan(ws))
              { w = wd; }
            else if (ws == 0)
              { w = 0; }
            else if (isnan(wd))
              { w = NAN; }
            else
              { w = ws*wd; }
            (*wdp) = (float)w;
          }
      }
     
    if (debug) { fprintf(stderr, "exiting %s ...\n", __FUNCTION__); }
  }
  
