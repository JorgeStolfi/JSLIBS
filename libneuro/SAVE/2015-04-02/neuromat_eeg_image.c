/* See {neuromat_{eeg_|}image.h}. */
/* Last edited on 2013-12-06 02:48:05 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <r2.h>
#include <float_image.h>
#include <float_pnm_image.h>
#include <float_image_paint.h>
#include <frgb.h>
#include <frgb_path.h>
#include <affirm.h>

#include <jsfile.h>

#include <neuromat_eeg.h>
#include <neuromat_image.h>
#include <neuromat_eeg_image.h>

float_image_t *neuromat_eeg_image_schematic_head_mask(int NX, int NY, r2_t *ctr, r2_t *rad)
  {
    float_image_t *img = float_image_new(1, NX, NY);
    float_image_fill(img, 0.0);
    float vfill = 1.0; /* Ellipse fill color. */
    float vdraw = NAN; /* Outline color. */
    float hwd = 0.0; /* Half-width of drawing pen. */
    int msub = 3; /* Subsampling order. */
    float_image_paint_ellipse_aligned(img, 0, ctr, rad, hwd, vfill, vdraw, msub);
    return img;
  }
  
void neuromat_eeg_image_draw_electrodes
  ( float_image_t *img,
    int ne, 
    r2_t pos[], 
    double drad, 
    int ilo, 
    float vlo[], 
    float vhi[]
  )
  { 
    int NC = (int)img->sz[0];
    int ie;
    for (ie = 0; ie < ne; ie++) 
      { r2_t *posi = &(pos[ie]);
        double cx = posi->c[0];
        double cy = posi->c[1];
        float *vfill = (ie == ilo ? vlo : vhi); /* Fill color vector. */
        float vdraw = NAN;   /* Outline color. (Not used) */
        float hwd = 0.0;     /* Half-width of drawing pen. */
        bool_t round = TRUE; /* Round dots? */
        bool_t diagonal = FALSE; /* Diagonal crosses? (Not used.) */
        int msub = 3; /* Subsampling order. */
        int c;
        for (c = 0; c < NC; c++)
          { float_image_paint_dot(img, c, cx, cy, drad, hwd, round, diagonal, vfill[c], vdraw, msub); } 
      }
  }

void neuromat_eeg_image_paint_potentials
  ( int ne, 
    double val[], 
    float_image_t *bas[],
    float_image_t *msk, 
    int c, float_image_t *img
  )
  {
    int NC = (int)img->sz[0];
    int NX = (int)img->sz[1];
    int NY = (int)img->sz[2];
    
    assert(msk->sz[0] == 1);
    assert(msk->sz[1] == NX);
    assert(msk->sz[2] == NY);
    
    demand((c >= 0) && (c < NC), "invalid channel");

    int ie;
    int ix, iy;
    for (iy = 0; iy < NY; iy++)
      { for (ix = 0; ix < NX; ix++)
          { double wxy = float_image_get_sample(msk, 0, ix, iy);
            if (wxy != 0)
              { double sum = 0;
                for (ie = 0; ie < ne; ie++) 
                  { double basi = float_image_get_sample(bas[ie], 0, ix, iy);
                    sum += val[ie] * basi;
                  }
                assert(! isnan(sum));
                float *smp = float_image_get_sample_address(img, c, ix, iy);
                double val = wxy*sum + (1-wxy)*(*smp);
                (*smp) = (float)val;
              }
          }
      }
  }

void neuromat_eeg_image_paint_marker_tics
  ( float_image_t *img, 
    int bxmin, int bxmax, 
    int mymin, int mymax, 
    int nt,
    int nc, 
    double **val, 
    int ic_mark
  )    
  {
    auto void paint_marker(double xa, double xb);
      /* Paint a marker from fractional column indices {xa} to {xb}. */
    
    /* Plot frames where marker is nonzero, merging any markers that are less than 1 pixel apart: */
    double xup = -INF; /* Fractional abscissa of last start-of-mark. */
    double xdn = -INF; /* Fractional abscissa of last end-of-mark. */
    int it;
    for (it = 0; it < nt; it++)
      { double vmark = val[it][ic_mark];
        demand(! isnan(vmark), "marker channel is NAN");
        if (vmark != 0)
          { double xlo = neuromat_image_slider_bar_position(bxmin, bxmax, (double)it, nt);
            double xhi = neuromat_image_slider_bar_position(bxmin, bxmax, (double)it + 1, nt);
            if (xhi - xlo < 1.0)
              { /* Fatten marker to 1 pixel: */
                double xmd = (xlo + xhi)/2;
                xlo = xmd - 0.5;
                xhi = xmd + 0.5;
              }
            if (xlo >= xdn + 1.0) 
              { /* Paint accumulated previous marker, if any: */
                if (xdn > xup) { paint_marker(xup, xdn); }
                /* Start new marker: */
                xup = xlo;
              }
            /* Current marker ends here so far: */
            xdn = xhi;
          }
      }
    /* Paint last marker, if any: */
    if (xdn > xup) { paint_marker(xup, xdn); }
    
    void paint_marker(double xa, double xb)
      {
        int NC = (int)img->sz[0];
        double hwd = 0;
        int c;
        for (c = 0; c < NC; c++) 
          { float vmark = 0.80f;
            float_image_paint_rectangle(img, c, xa, xb, (double)mymin, (double)mymax+1, hwd, vmark, NAN, 6);
          }
      }
  }

