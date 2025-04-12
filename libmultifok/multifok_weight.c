/* See {multifok_weight.h}. */
/* Last edited on 2025-04-11 09:05:42 by stolfi */

#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <r3.h>
#include <float_image.h>

#include <multifok_weight.h>

float_image_t *multifok_weight_from_height_and_normal(float_image_t *htd, double dhMax, float_image_t *nrm, double dnMax)
  {
    demand((htd != NULL) || (nrm != NULL), "at least one of {htd} and {nrm} must be given");
    
    int32_t NC_h = -1, NC_n = -1, NX = -1, NY = -1;
    if (htd != NULL)
      { float_image_get_size(htd, &NC_h, &NX, &NY);
        if (nrm != NULL)
          { NC_n = (int32_t)(nrm->sz[0]);
            float_image_check_size(nrm, -1, NX, NY, "inconsistent sizes {htd,nrm}");
          }
      }
    else
      { float_image_get_size(nrm, &NC_h, &NX, &NY); }
      
    if (htd != NULL) { demand(NC_h >= 1, "invalid {htd} channel count"); }
    if (nrm != NULL) { demand(NC_n >= 3, "invalid {nrm} channel count"); }
        
    auto double weight_from_htd(int32_t x, int32_t y);
      /* Return the weight for pixel {[x,y]} (either 0 or 1) considering {htd} alone.
        In particular, if {htd} is null, returns 1. */
        
    auto double weight_from_nrm(int32_t x, int32_t y);
      /* Return the weight for pixel {[x,y]} (either 0 or 1) considering {nrm} alone.
        In particular, if {nrm} is null, returns 1. */
        
    auto bool_t get_normal(int32_t x, int32_t y, r3_t *u_P);
      /* Tries to gets the normal vector {u=nrm[0..2,x,y]} and normalize
        it to unit length. If the pixel {[x,y]} does not exist or any
        coordinate of {u} is not finite or the norm of {u} is too close
        to zero, returns {FALSE} and {*u_P} is undefined. Othwerwise
        stores the normalized {u} into {*u_P} and returns {TRUE} */

    float_image_t *wht = float_image_new(1, NX, NY);
    for (int32_t y = 0; y < NY;y++)
      { for (int32_t x = 0; x < NX; x++)
          { double w = weight_from_htd(x, y);
            if (w != 0) { w *= weight_from_nrm(x, y); }
            float_image_set_sample(wht, 0, x, y, (float)w);
          }
      }
      
    return wht;
    
    double weight_from_htd(int32_t x, int32_t y)
      { if (htd == NULL) { return 1.0; }
        
        /* Check height deviation: */
        float dh = float_image_get_sample(htd, 0, x, y);
        if (! isfinite(dh)) { return 0.0; }
        demand(dh >= 0, "deviation cannot be negative");
        if (dh > dhMax) { return 0.0; }
        
        if (NC_h >= 2)
          { /* Get the weight from {htd} too: */
            float wh = float_image_get_sample(htd, 1, x, y);
            demand(isfinite(wh) && (wh >= 0), "invalid weight value in height deviation map");
            if (wh == 0) { return 0.0; }
          }
        return 1.0;
      }
    
    double weight_from_nrm(int32_t x, int32_t y)
      { if (nrm == NULL) { return 1.0; }
        /* Get normal at pixel: */
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
        if (NC_n >= 4)
          { /* Get the weight from {nrm} too: */
            float wn = float_image_get_sample(nrm, 1, x, y);
            demand(isfinite(wn) && (wn >= 0), "invalid weight value in normal map");
            if (wn == 0) { return 0.0; }
          }
        return 1.0;
      }
        
    bool_t get_normal( int32_t x, int32_t y, r3_t *u_P)
      { if ((x < 0) || (x >= NX) || (y < 0) || (y >= NY)) { return FALSE; }
        r3_t u;
        for (int32_t c = 0; c < 3; c++) 
          { u.c[c] = float_image_get_sample(nrm, c, x, y);
            if (! isfinite(u.c[c])) { return FALSE; } 
          }
        /* Try to normalize, just in case: */
        double u_mag = r3_dir(&u, &u);
        if (u_mag < 1.0e-4) { return FALSE; }
        (*u_P) = u;
        return TRUE;
      }
  }
