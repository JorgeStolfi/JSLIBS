/* See pst_shading.h */
/* Last edited on 2024-12-22 22:52:27 by stolfi */

#include <stdio.h>
#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <float_image.h>
#include <r3.h> 
#include <affirm.h>
#include <argparser.h>

#include <pst_basic.h>
#include <pst_lamp.h>
#include <pst_light.h>
#include <pst_normal_map.h>
#include <pst_shading.h>

float_image_t *pst_shading_difference_image
  ( float_image_t *NRM, 
    float_image_t *AIMG, 
    float_image_t *BIMG
  )
  { int32_t NC, NX, NY;
    float_image_get_size(AIMG, &NC, &NX, &NY);
    float_image_check_size(BIMG, NC, NX, NY);
    float_image_check_size(NRM, 3, NX, NY);

    float_image_t *DIF = float_image_new(NC, NX, NY);
    
    int32_t c, x, y;
    for (y = 0; y < NY; y++)
      { for (x = 0; x < NX; x++)
          { r3_t nrm = pst_normal_map_get_pixel(NRM, x, y);
            bool_t ok = (r3_L_inf_norm(&nrm) != 0); 
            for (c = 0; c < NC; c++) 
              { float sd;
                if (ok) 
                  { float sa = float_image_get_sample(AIMG, c, x, y); 
                    float sb = float_image_get_sample(BIMG, c, x, y); 
                    sd = sa - sb;
                  }
                else
                  { sd = 0.0; }
                float_image_set_sample(DIF, c, x, y, sd); 
              }
          }
      }
    return DIF;
  }    

void pst_shading_add_diffuse_single
  ( float_image_t *NRM, 
    float_image_t *CLR,
    pst_lamp_t *src,
    float_image_t *IMG
  )
  { /* Create a vector with that single lamp: */
    pst_lamp_vec_t svec = pst_lamp_vec_new(1);
    svec.e[0] = src;
    pst_light_t lht = pst_light_from_lamps(svec);
    /* Call the general Lambertian shading procedure with it: */
    pst_shading_add_diffuse(NRM, CLR, &lht, IMG);
    /* Free the vector: */
    free(svec.e);
  }

void pst_shading_add_diffuse
  ( float_image_t *NRM, 
    float_image_t *CLR,
    pst_light_t *lht,
    float_image_t *IMG
  )
  { pst_lamp_vec_t *lmpv = &(lht->lmpv);
    uint32_t NS = lmpv->ne;

    /* Get and check image sizes: */
    int32_t NC, NX, NY;
    float_image_get_size(IMG, &NC, &NX, &NY);
    float_image_check_size(NRM, 3, NX, NY);
    if (CLR != NULL) { float_image_check_size(CLR, NC, NX, NY); }
    
    double clr[NC]; /* Intrinsic color of scene at some pixel. */
    double lum[NC]; /* Diffuse illumination due to all lamps. */

    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { r3_t nrm = pst_normal_map_get_pixel(NRM, x, y);
            if (r3_L_inf_norm(&nrm) != 0)
              { /* Valid pixel, compute its shaded color. */
                /* Clear apparent color, get intrinsic color, check if black: */
                bool_t black = TRUE;
                for (uint32_t c = 0; c < NC; c++) 
                  { lum[c] = 0.0;
                    clr[c] = (CLR == NULL ? 1.0 : float_image_get_sample(CLR, (int32_t)c, x, y)); 
                    if (clr[c] != 0.0) { black = FALSE; }
                  }
                if (!black)
                  { /* Shine all lamps on it: */
                    for (uint32_t i = 0; i < NS; i++)
                      { pst_lamp_t *src = lmpv->e[i];
                        r3_t *dir = &(src->dir);
                        double crad = src->crad;
                        /* Compute its geometric factor for diffusion term: */
                        double coef = pst_lamp_geom_factor(&nrm, dir, crad);
                        /* Accumulate light on each channel: */
                        double_vec_t *pwr = &(src->pwr);
                        assert(pwr->ne == NC); 
                        for (uint32_t c = 0; c < NC; c++) { lum[c] += coef*pwr->e[c]; }
                      }
                    /* Add color to {IMG}: */
                    for (uint32_t c = 0; c < NC; c++) 
                      { float *p = float_image_get_sample_address(IMG, (int32_t)c, x, y); 
                        (*p) += (float)(lum[c]*clr[c]);
                      }
                  }
              }
          }
      }
  }

