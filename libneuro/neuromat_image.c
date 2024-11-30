/* See {neuromat_image.h}. */
/* Last edited on 2021-08-29 00:21:32 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <r2.h>
#include <sign.h>
#include <float_image.h>
#include <float_image_paint.h>
#include <frgb.h>
#include <frgb_ops.h>
#include <frgb_path.h>
#include <affirm.h>

#include <jsfile.h>

#include <neuromat_eeg.h>
#include <neuromat_image.h>

void neuromat_image_colorize_field
  ( float_image_t *cim,
    float_image_t *fld, 
    float_image_t *msk, 
    double vmax, 
    int32_t style
  )
  { 
    demand(fld->sz[0] == 1, "image should be monochromatic");
    int32_t NX = (int32_t)fld->sz[1];
    int32_t NY = (int32_t)fld->sz[2];
    
    if (msk != NULL)
      { demand(msk->sz[0] == 1,  "mask must be monochromatic");
        demand(msk->sz[1] == NX, "mask with wrong {NX}");
        demand(msk->sz[2] == NY, "mask with wrong {NY}");
      }
    
    int32_t NC = (int32_t)cim->sz[0];
    demand(NC == 4, "result image must be RGBA");
    demand(cim->sz[1] == NX, "result with wrong {NX}");
    demand(cim->sz[2] == NY, "result with wrong {NY}");
    
    demand(vmax > 0, "invalid {vmax}");

    for (uint32_t iy = 0;  iy < NY; iy++)
      { for (uint32_t ix = 0;  ix < NX; ix++)
          { /* Get mask and image value: */
            double fxy = (double)float_image_get_sample(fld, 0, ix, iy); 
            float mxy = (msk == NULL ? 1.000f : float_image_get_sample(msk, 0, ix, iy)); 
            demand(! isnan(mxy), "mask is NAN");
            if ((mxy <= 0) || (isnan(fxy)))
              { /* Set pixel to {bgr} with opacity zero: */
                float_image_fill_pixel(cim, ix, iy, 0.000f);
              }
            else
              { /* Convert {fld} from {[-vmax _ +vmax]} to {[-1 _ +1]}: */
                double z = fxy/vmax;
                if (z < -1) { z = -1; }
                if (z > +1) { z = +1; }
                frgb_t rgb = frgb_path_map_signed(z, 1, style);
                /* Set pixel to {rgb} with opacity {mxy}: */
                float ccolor[NC];  for (uint32_t c = 0;  c < NC; c++) { ccolor[c] = (c < 3 ? rgb.c[c] : mxy); }
                /* Save colorized pixel: */
                float_image_set_pixel(cim, ix, iy, ccolor);
              }
          }
      }
  }

void neuromat_image_colorize_signed_overlay
  ( float_image_t *cim,
    float_image_t *ovr, 
    sign_t sgn,
    frgb_t fc
  )
  {
    int32_t NC = (int32_t)cim->sz[0];
    demand(NC == 4, "result image must be RGBA");
    int32_t NX = (int32_t)cim->sz[1];
    int32_t NY = (int32_t)cim->sz[2];
    
    demand(ovr->sz[0] == 1,  "overlay must be monochromatic");
    demand(ovr->sz[1] == NX, "overlay with wrong {NX}");
    demand(ovr->sz[2] == NY, "overlay with wrong {NY}");

    /* Add opacity to {fc}: */
    float color[NC]; for (uint32_t c = 0;  c < NC; c++) { color[c] = (c < 3 ? fc.c[c] : 1.000f); }
    
    for (uint32_t iy = 0;  iy < NY; iy++)
      { for (uint32_t ix = 0;  ix < NX; ix++)
          { /* Get overlay value {oxy}: */
            float oxy = float_image_get_sample(ovr, 0, ix, iy); 
            if (isnan(oxy) || (oxy*((float)sgn) <= 0))
              { /* Set the pixel to transparent black: */
                float_image_fill_pixel(cim, ix, iy, 0.000f);
              }
            else 
              { /* Set the pixel to semitransparent {fc}: */
                color[3] = fabsf(oxy);
                float_image_set_pixel(cim, ix, iy, color);
              }
          }
      }
  }

void neuromat_image_paint_time_track
  ( float_image_t *img,
    int32_t hw, 
    int32_t xlo, 
    int32_t xsz, 
    int32_t y, 
    frgb_t fc
  )
  {
    int32_t NC = (int32_t)img->sz[0];
    int32_t NX = (int32_t)img->sz[1];
    int32_t NY = (int32_t)img->sz[2];
    demand(NC == 4, "result image should be RGBA");
    demand((xlo >= 0) && (xlo + xsz <= NX), "invalid timeline abscissas");
    demand((y >= hw) && (y <= NY - hw), "invalid timeline ordinate");
    
    int32_t xmin = xlo;
    int32_t xmax = xlo + xsz - 1;
    int32_t ymin = y - hw;
    int32_t ymax = y + hw - 1;
    /* Add opacity to {fc}: */
    float color[NC]; for (uint32_t c = 0;  c < NC; c++) { color[c] = (c < 3 ? fc.c[c] : 1.000f); }
    float_image_fill_rectangle_pixels(img, xmin, xmax, ymin, ymax, color);
  }

void neuromat_image_paint_time_range_and_tics
  ( float_image_t *img,
    int32_t hw, 
    double tlo,
    double tini,
    double tfin,
    double thi, 
    int32_t xlo, 
    int32_t xsz,
    int32_t y, 
    frgb_t fc
  )
  {
    int32_t NC = (int32_t)img->sz[0];
    demand(NC == 4, "result image should be RGBA");
    if (tini > tfin)
      { /* Empty time interval: */
        return;
      }
    else if (tini == tfin)
      { /* Single instant: */
        neuromat_image_paint_tic(img, hw, tlo, tini, thi, xlo, xsz, y, fc);
      }
    else
      { /* Non-empty interval: */
        neuromat_image_paint_time_range(img, hw, tlo, tini, tfin, thi, xlo, xsz, y, fc);
        neuromat_image_paint_tic(img, hw, tlo, tini, thi, xlo, xsz, y, fc);
        neuromat_image_paint_tic(img, hw, tlo, tfin, thi, xlo, xsz, y, fc);
      }
  }

void neuromat_image_paint_time_range
  ( float_image_t *img,
    int32_t hw, 
    double tlo,
    double tini,
    double tfin,
    double thi, 
    int32_t xlo, 
    int32_t xsz,
    int32_t y, 
    frgb_t fc
  )
  { 
    int32_t NC = (int32_t)img->sz[0];
    demand(NC == 4, "result image should be RGBA");
    /* Provide opacity channel 3 for {fc} if needed: */
    float vfill[NC]; for (uint32_t c = 0;  c < NC; c++) { vfill[c] = (c < 3 ? fc.c[c] : 1.000f); }
    
    demand(thi - tlo >= 1.0e-10, "global interval {tlo_thi} too small");

    if (tini < tlo) { tini = tlo; }
    if (tini > thi) { tini = thi; }
    
    double scale = ((double)xsz)/(thi - tlo);
    int32_t xmin = xlo + (int32_t)floor(scale*(tini - tlo) + 0.5);
    int32_t xmax = xlo + (int32_t)floor(scale*(tfin - tlo) + 0.5) - 1;
    if (xmax < xmin) { return; }
    int32_t ymin = y - hw;
    int32_t ymax = y + hw - 1;
    float_image_fill_rectangle_pixels(img, xmin, xmax, ymin, ymax, vfill); 
  }

void neuromat_image_paint_tic
  ( float_image_t *img,
    int32_t hw, 
    double tlo,
    double t,
    double thi, 
    int32_t xlo, 
    int32_t xsz,
    int32_t y, 
    frgb_t fc
  )
  {
    int32_t NC = (int32_t)img->sz[0];
    demand(NC == 4, "result image should be RGBA");
    /* Provide opacity channel 3 for {fc} if needed: */
    float vfill[NC]; for (uint32_t c = 0;  c < NC; c++) { vfill[c] = (c < 3 ? fc.c[c] : 1.000f); }
    
    demand(thi - tlo >= 1.0e-10, "global interval {tlo_thi} too small");

    if ((t < tlo) || (t > thi)) { return; }
    
    double scale = ((double)xsz)/(thi - tlo);
    int32_t xctr = xlo + (int32_t)floor(scale*(t - tlo) + 0.5);
    int32_t xmin = xctr - hw;
    int32_t xmax = xctr + hw - 1;
    int32_t ymin = y - 2*hw;
    int32_t ymax = y + 2*hw - 1;
    float_image_fill_rectangle_pixels(img, xmin, xmax, ymin, ymax, vfill);
  }
  
void neuromat_image_paint_slider
  ( float_image_t *img, 
    int32_t hw,
    int32_t hh, 
    double tlo,
    double t,
    double thi, 
    int32_t xlo, 
    int32_t xsz,
    int32_t ylo,
    int32_t yhi, 
    frgb_t fc
  )
  { 
    int32_t NC = (int32_t)img->sz[0];
    demand(NC == 4, "result image should be RGBA");
    /* Provide opacity channel 3 for {fc} if needed: */
    float vfill[NC]; for (uint32_t c = 0;  c < NC; c++) { vfill[c] = (c < 3 ? fc.c[c] : 1.000f); }
    
    demand(thi - tlo >= 1.0e-10, "global interval {tlo_thi} too small");

    if ((t < tlo) || (t > thi)) { return; }
    
    double scale = ((double)xsz)/(thi - tlo);
    int32_t xctr = xlo + (int32_t)floor(scale*(t - tlo) + 0.5);

    int32_t xmin = xctr - hw;
    int32_t xmax = xctr + hw - 1;
    
    /* Paint the lower slider: */
    int32_t ymin = ylo - hh;
    int32_t ymax = ylo - 1;
    float_image_fill_rectangle_pixels(img, xmin, xmax, ymin, ymax, vfill);
    
    /* Paint the upper slider: */
    ymin = yhi;
    ymax = yhi + hh - 1;
    float_image_fill_rectangle_pixels(img, xmin, xmax, ymin, ymax, vfill);
  }
