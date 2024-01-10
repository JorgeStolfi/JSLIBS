/* See {neuromat_image.h}. */
/* Last edited on 2013-12-06 04:50:26 by stolfilocal */

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

float_image_t *neuromat_image_make_slider_bar(int NC, int NX, int NY, int *bxminP, int *bxmaxP)    
  {
    demand(NY % 2 == 0, "slider bar image height must be multiple of 2");
    demand(NY >= 4, "slider bar image height too small");
    demand(NX > NY, "slider bar image too narrow");
    
    /* Size and position of trigger axis bar (pixels): */
    int bysize = 2;
    int bymin = (NY - bysize)/2;
    int bymax = NY - 1 - bymin;
    
    /* Paint the axis bar: */
    int bxmin = bymin;
    int bxmax = NX - 1 - bxmin;

    float_image_t *img = float_image_new(NC, NX, NY);
    float_image_fill(img, 0.0);
    float vbar = 0.50;
    float_image_fill_rectangle(img, bxmin, bxmax, bymin, bymax, vbar);
    (*bxminP) = bxmin;
    (*bxmaxP) = bxmax;
    return img;
  }
    
double neuromat_image_slider_bar_position(int bxmin, int bxmax, double t, int nt)
  { 
    /* Linear interpolation along slider bar: */
    demand(bxmin <= bxmax, "bad image column range");
    double bxsize = bxmax + 1 - bxmin;
    return bxmin + (t*bxsize)/nt;
  }

void neuromat_image_paint_slider_dot(float_image_t *img, int bxmin, int bxmax, double yctr, double t, int nt)
  { double xctr = neuromat_image_slider_bar_position(bxmin, bxmax, t, nt);
    double rad = 3.0;
    double hwd = 1.0;
    bool_t round = TRUE;
    bool_t diagonal = FALSE;
    int NC = (int)img->sz[0];
    int c;
    for (c = 0; c < NC; c++) 
      { float vfill = NAN;
        float vdraw = (float)(NC == 1 ? 1.0 : 1.0 - c/(NC - 0.9));
        float_image_paint_dot(img, c, xctr, yctr, rad, hwd, round, diagonal, vfill, vdraw, 6);
      }
  }

float_image_t *neuromat_image_colorize(float_image_t *fim, float_image_t *msk, double vmax, int style, frgb_t *bgr)
  { 
    demand(fim->sz[0] == 1, "image should be monochromatic");
    int NX = (int)fim->sz[1];
    int NY = (int)fim->sz[2];
    
    if (msk != NULL)
      { demand(msk->sz[0] == 1, "mask must be monochromatic");
        demand(msk->sz[1] == NX, "mask with wrong {NX}");
        demand(msk->sz[2] == NY, "mask with wrong {NY}");
      }
    
    float_image_t *cim = float_image_new(3, NX, NY);

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

    return cim;
  }

void neuromat_image_paint_overlay(float_image_t *fim, float_image_t *ovr, float col[])
  {
    int NC = (int)fim->sz[0];
    int NX = (int)fim->sz[1];
    int NY = (int)fim->sz[2];
    
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
                    float_image_set_pixel(fim, ix, iy, col);
                  }
                else
                  { /* Overlay is semi-transparent, mix colors: */
                    float fxy[NC]; 
                    float_image_get_pixel(fim, ix, iy, fxy);
                    for (ic = 0; ic < NC; ic++) { fxy[ic] = (float)(oxy*col[ic] + (1-oxy)*fxy[ic]); }
                    float_image_set_pixel(fim, ix, iy, fxy);
                  }
              }
          }
      }
  }
