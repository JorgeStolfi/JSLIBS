/* See {neuromat_image.h}. */
/* Last edited on 2020-10-11 01:46:31 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <r2.h>
#include <float_image.h>
#include <float_image_paint.h>
#include <frgb.h>
#include <frgb_ops.h>
#include <frgb_path.h>
#include <affirm.h>

#include <jsfile.h>

#include <neuromat_eeg.h>
#include <neuromat_image.h>

void neuromat_image_paint_time_track
  ( float_image_t *img, 
    int hw, 
    int xlo, 
    int xsz, 
    int y, 
    frgb_t *fc
  )
  {
    int NC = (int)img->sz[0];
    int NX = (int)img->sz[1];
    int NY = (int)img->sz[2];
    demand(NC == 3, "image should be colored");
    demand((xlo >= 0) && (xlo + xsz <= NX), "invalid timeline abscissas");
    demand((y >= hw) && (y <= NY - hw), "invalid timeline ordinate");
    
    int xmin = xlo;
    int xmax = xlo + xsz - 1;
    int ymin = y - hw;
    int ymax = y + hw - 1;
    float_image_fill_rectangle_pixels(img, xmin, xmax, ymin, ymax, fc->c);
  }

void neuromat_image_paint_marker_ranges
  ( float_image_t *img, 
    int hw, 
    int nt, 
    int nc,
    double **val,
    int ic_mark, 
    int xlo, 
    int xsz,
    int y,
    frgb_t *fc
  )
  { double tlo = 0.0; /* Time just before first frame. */
    double thi = (double)nt; /* Time just after last frame. */
    /* Scan samples of channel {ic}: */
    /* Pretend that marker channel is zero before and after all frames. */
    int it_ini = -1;  /* Index of start frame of run, or -1 if not started yet. */
    int it_fin = -1;  /* Index of end frame of run, or -1 if not started yet. */
    for (int it = 0; it <= nt; it++)
      { double smp = (it >= nt ? 0.0 : val[it][ic_mark]);
        assert(! isnan(smp));
        if (smp != 0.0)
          { /* Start or extend a run: */
            if (it_ini < 0) { it_ini = it; }
            it_fin = it;
          }
        else if (it_ini >= 0)
          { /* End of run: */
            assert((it_ini >= 0) && (it_ini <= it_fin) && (it_fin < nt));
            double tini = it_ini + 0.5;
            double tfin = it_fin + 0.5;
            neuromat_image_paint_time_range_and_tics(img, hw, tlo, tini, tfin, thi, xlo, xsz, y, fc);
            it_ini = -1;
            it_fin = -1;
          }
      }
  }

void neuromat_image_paint_time_range_and_tics
  ( float_image_t *img, 
    int hw, 
    double tlo,
    double tini,
    double tfin,
    double thi, 
    int xlo, 
    int xsz,
    int y,
    frgb_t *fc
  )
  {
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
    int hw, 
    double tlo,
    double tini,
    double tfin,
    double thi, 
    int xlo, 
    int xsz,
    int y,
    frgb_t *fc
  )
  { 
    int NC = (int)img->sz[0];
    demand(NC == 3, "image should be colored");
    
    demand(thi - tlo >= 1.0e-10, "global interval {tlo_thi} too small");

    if (tini < tlo) { tini = tlo; }
    if (tini > thi) { tini = thi; }
    
    double scale = ((double)xsz)/(thi - tlo);
    double xmin = xlo + scale*(tini - tlo);
    double xmax = xlo + scale*(tfin - tlo);
    if (xmax - xmin < 1.0e-2) { return; }
    double ymin = (double)(y - hw);
    double ymax = (double)(y + hw);
    float vdraw = NAN;
    for (int c = 0; c < NC; c++)
      { float vfill = fc->c[c];
        (void)float_image_paint_rectangle(img, c, xmin, xmax, ymin, ymax, 0.0, vfill, vdraw, 6);
      }
  }

void neuromat_image_paint_tic
  ( float_image_t *img, 
    int hw, 
    double tlo,
    double t,
    double thi, 
    int xlo, 
    int xsz,
    int y,
    frgb_t *fc
  )
  {
    int NC = (int)img->sz[0];
    demand(NC == 3, "image should be colored");
    
    demand(thi - tlo >= 1.0e-10, "global interval {tlo_thi} too small");

    if ((t < tlo) || (t > thi)) { return; }
    
    double scale = ((double)xsz)/(thi - tlo);
    double xctr = xlo + scale*(t - tlo);
    double xmin = xctr - (double)hw;
    double xmax = xctr + (double)hw;
    double ymin = (double)(y - 2*hw);
    double ymax = (double)(y + 2*hw);
    float vdraw = NAN;
    for (int c = 0; c < NC; c++)
      { float vfill = fc->c[c];
        (void)float_image_paint_rectangle(img, c, xmin, xmax, ymin, ymax, 0.0, vfill, vdraw, 6);
      }
  }
  
void neuromat_image_paint_slider
  ( float_image_t *img, 
    int hw,
    double tlo,
    double t,
    double thi, 
    int xlo, 
    int xsz,
    int ylo,
    int yhi,
    frgb_t *fc
  )
  { 
    int NC = (int)img->sz[0];
    demand(NC == 3, "image should be colored");
    
    demand(thi - tlo >= 1.0e-10, "global interval {tlo_thi} too small");

    if ((t < tlo) || (t > thi)) { return; }
    
    double scale = ((double)xsz)/(thi - tlo);
    double xctr = xlo + scale*(t - tlo);

    double xmin = xctr - (double)hw;
    double xmax = xctr + (double)hw;
    
    for (int ik = 0; ik < 2; ik++)
      { double yctr = (double)(ik == 0 ? ylo - 4*hw : yhi + 4*hw);
        double ymin = yctr - 2*hw;
        double ymax = yctr + 2*hw;
        float vdraw = NAN;
        for (int c = 0; c < NC; c++)
          { float vfill = fc->c[c];
            (void)float_image_paint_rectangle(img, c, xmin, xmax, ymin, ymax, 0.0, vfill, vdraw, 6);
          }
      }
  }

void neuromat_image_colorize
  ( float_image_t *fim, 
    float_image_t *msk, 
    double vmax, 
    int style, 
    frgb_t *bgr,
    float_image_t *cim
  )
  { 
    demand(fim->sz[0] == 1, "image should be monochromatic");
    int NX = (int)fim->sz[1];
    int NY = (int)fim->sz[2];
    
    if (msk != NULL)
      { demand(msk->sz[0] == 1,  "mask must be monochromatic");
        demand(msk->sz[1] == NX, "mask with wrong {NX}");
        demand(msk->sz[2] == NY, "mask with wrong {NY}");
      }
    
    demand(cim->sz[0] == 3,  "result must be RGB");
    demand(cim->sz[1] == NX, "result with wrong {NX}");
    demand(cim->sz[2] == NY, "result with wrong {NY}");

    int ix, iy;
    for (iy = 0; iy < NY; iy++)
      { for (ix = 0; ix < NX; ix++)
          { /* Get mask and image value: */
            double fxy = (double)float_image_get_sample(fim, 0, ix, iy); 
            double mxy = (msk == NULL ? 1.0 : (double)float_image_get_sample(msk, 0, ix, iy)); 
            if ((mxy <= 0) || (isnan(fxy)))
              { float_image_set_pixel(cim, ix, iy, bgr->c); }
            else
              { /* Remove mask factor: */
                fxy /= mxy;
                double z = fxy/vmax;
                if (z < -1) { z = -1; }
                if (z > +1) { z = +1; }
                frgb_t rgb = frgb_path_map_signed(z, 1, style);
                /* Put back mask factor assuming background is zero: */
                frgb_mix(mxy, &rgb, 1-mxy, bgr);
                /* Save colorized pixel: */
                float_image_set_pixel(cim, ix, iy, rgb.c);
              }
          }
      }
  }

void neuromat_image_paint_overlay(float_image_t *img, float_image_t *ovr, float col[])
  {
    int NC = (int)img->sz[0];
    int NX = (int)img->sz[1];
    int NY = (int)img->sz[2];
    
    demand(ovr->sz[0] == 1,  "overlay must be monochromatic");
    demand(ovr->sz[1] == NX, "overlay with wrong {NX}");
    demand(ovr->sz[2] == NY, "overlay with wrong {NY}");

    int ic, ix, iy;
    for (iy = 0; iy < NY; iy++)
      { for (ix = 0; ix < NX; ix++)
          { /* Get overlay opacity {oxy}: */
            double oxy = (double)float_image_get_sample(ovr, 0, ix, iy); 
            demand(! isnan(oxy), "overlay is NAN");
            if (oxy != 0)
              { demand((oxy >= 0) && (oxy <= 1.0), "overlay sample out of range");
                if (oxy >= 1)
                  { /* Overlay is opaque: */
                    float_image_set_pixel(img, ix, iy, col);
                  }
                else
                  { /* Overlay is semi-transparent, mix colors: */
                    float fxy[NC]; 
                    float_image_get_pixel(img, ix, iy, fxy);
                    for (ic = 0; ic < NC; ic++) { fxy[ic] = (float)(oxy*col[ic] + (1-oxy)*fxy[ic]); }
                    float_image_set_pixel(img, ix, iy, fxy);
                  }
              }
          }
      }
  }
