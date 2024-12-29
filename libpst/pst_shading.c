/* See pst_shading.h */
/* Last edited on 2024-12-29 02:06:44 by stolfi */

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

void pst_shading_pixel_color_and_alpha
  ( int32_t x,
    int32_t y,
    pst_normal_func_t *normal,
    pst_albedo_func_t *albedo,
    pst_shading_func_t* shading,
    frgb_t *val_P,
    float *alpha_P
  );
  /* Computes the mean apparent color (per-channel radiance) {val} and
    opacity {alpha} of the synthetic image in the pixel with bottom left
    corner {(x,y)}, as a weighted average of the color {val(p,c)} in
    channel {c} at several points {p} of the image in and around the
    pixel.
    
    The point color {val} takes into account the surface normal
    direction {normal(p)}, the surface albedo {albedo(p)}, and the
    lighting of the surface at point {p}, embodied in the function
    {shading(nrm)}.
    
    A point is omitted from the average if any of these functions returns
    {(NAN,NAN,NAN)}.  The opacity {alpha} will be the fraction of
    total weight that got included in the average. */

/* IMPLEMENTATIONS */

float_image_t* pst_shading_make_image
  ( int32_t NC, int32_t NX, int32_t NY,
    pst_normal_func_t *normal,
    pst_albedo_func_t *albedo,
    pst_shading_func_t* shading
  )
  { demand((NC == 3) || (NC == 4), "invalid channel count {NC}");
    float_image_t *img = float_image_new(NC, NX, NY);
    pst_shading_paint(img, 0, NX, 0, NY, normal, albedo, shading);
    return img;
  }

double pst_shading_paint
  ( float_image_t *img,
    int32_t xlo, int32_t xhi,
    int32_t ylo, int32_t yhi,
    pst_normal_func_t *normal,
    pst_albedo_func_t *albedo,
    pst_shading_func_t* shading
  )
  { int32_t NC, NX, NY;
    float_image_get_size(img, &NC, &NX, &NY);
    demand((NC == 3) || (NC == 4), "image should have 3 or 4 channels");

    if (xlo < 0) { xlo = 0; }
    if (xhi >= NX) { xhi = NX-1; }

    if (ylo < 0) { ylo = 0; }
    if (yhi >= NY) { yhi = NY-1; }
    
    double sum_alpha = 0.0;
    
    if ((xlo > xhi) || (ylo > yhi)) { return sum_alpha; }

    for (int32_t x = xlo; x <= xhi; x++)
      { for (int32_t y = ylo; y <= yhi; y++)
          { frgb_t val;
            float alpha;
            pst_shading_pixel_color_and_alpha
              ( x, y, normal, albedo, shading, &val, &alpha);
            for (int32_t c = 0; c < NC; c++)
              { if (alpha > 0)
                  { float *smp = float_image_get_sample_address(img, c, x, y);
                    if (c < 3)
                      { (*smp) = (1-alpha)*(*smp) + val.c[c]; }
                    else
                      { assert((c == 3) && (NC == 4));
                        (*smp) = alpha;
                      }
                    sum_alpha += alpha;
                  }
              }
          } 
      }
    return sum_alpha;
  }

void pst_shading_pixel_color_and_alpha
  ( int32_t x,
    int32_t y,
    pst_normal_func_t *normal,
    pst_albedo_func_t *albedo,
    pst_shading_func_t* shading,
    frgb_t *val_P,
    float *alpha_P
  )
  { /* !!! Should use Hann window weights extending outside the pixel. !!! */
    int32_t HS = 2;  /* Pixel subsamples are numbered {-HS..+HS}. */
    double tot_w = (2*HS+1)*(2*HS+1); /* Total weright of all sampoints. */
    double sum_wv[3] = { 0,0,0 }; /* Sum of {wt*val} for valid pts. */
    double sum_w = 0; /* Sum of {wt} for valid pixs. */
    for (int32_t kx = -HS; kx <= +HS; kx++)
      { double px = (double)x + 0.5 + (double)(kx)/(HS+1.0);
        for (int32_t ky = -HS; ky <= +HS; ky++)
          { double py = (double)y + 0.5 + 0.5*(double)(ky)/(HS+1.0);
            double wt = 1.0;
            r2_t p = (r2_t){{ px, py }}; /* Subsampling point. */
            r3_t nrm = normal(&p);
            frgb_t alb = albedo(&p);
            frgb_t shd = shading(&nrm); 
            if (! (isnan(nrm.c[0]) || isnan(alb.c[0]) || isnan(shd.c[0])))
              { for (int32_t c = 0; c < 3; c++)
                  { sum_wv[c] += wt * shd.c[c]*alb.c[c]; }
                sum_w += wt;
              }
           }
       }
     (*alpha_P) = (float)(sum_w / tot_w);
     if (sum_w == 0)
       { (*val_P) = frgb_Black; }
     else
       { for (int32_t c = 0; c < 3; c++) 
           { val_P->c[c] = (float)(sum_wv[c] / sum_w); }
       }
  }

float_image_t *pst_shading_difference_image
  ( float_image_t *AIMG, 
    float_image_t *BIMG, 
    float_image_t *MSK, 
    float_image_t *NRM
  )
  { int32_t NC, NX, NY;
    float_image_get_size(AIMG, &NC, &NX, &NY);
    float_image_check_size(BIMG, NC, NX, NY);
    if (NRM != NULL) { float_image_check_size(NRM, 3, NX, NY); }
    if (MSK != NULL) { float_image_check_size(NRM, 1, NX, NY); }

    float_image_t *DIF = float_image_new(NC, NX, NY);
    
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { float msk = (MSK == NULL ? 1.0f : float_image_get_sample(MSK, 0, x, y));
            demand(isfinite(msk), "invalid mask value");
            r3_t nrm = pst_normal_map_get_pixel(NRM, x, y);
            double dn = r3_L_inf_norm(&nrm);
            bool_t ok = ((msk >= 0.5) && isfinite(dn) && (dn > 1.0e-16)); 
            for (int32_t c = 0; c < NC; c++) 
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
    pst_light_t *lht = pst_light_from_lamps(1, svec);
    /* Call the general Lambertian shading procedure with it: */
    pst_shading_add_diffuse(NRM, CLR, lht, IMG);
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
                for (int32_t c = 0; c < NC; c++) 
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
                        frgb_t *pwr = &(src->pwr);
                        for (int32_t c = 0; c < NC; c++) { lum[c] += coef*pwr->c[c]; }
                      }
                    /* Add color to {IMG}: */
                    for (int32_t c = 0; c < NC; c++) 
                      { float *p = float_image_get_sample_address(IMG, (int32_t)c, x, y); 
                        (*p) += (float)(lum[c]*clr[c]);
                      }
                  }
              }
          }
      }
  }

