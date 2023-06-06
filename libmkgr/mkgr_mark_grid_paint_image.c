/* See {mkgr_mark_grid_paint_image.h} */
/* Last edited on 2023-04-23 11:10:48 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <vec.h>
#include <r2.h>
#include <affirm.h>
#include <frgb.h>
#include <float_image_paint.h>

#include <mkgr_mark.h>
#include <mkgr_mark_grid.h>
#include <mkgr_mark_grid_paint_image.h>

void mkgr_mark_grid_paint_image
  ( float_image_t *img,
    int32_t chns,
    int32_t ch[],
    mkgr_mark_grid_t *gr, 
    double scale,
    r2_t *org, 
    int32_t m
  )
  {
    int32_t NC = (int32_t)img->sz[0];

    demand(chns >= 0, "invalid channel count");
    if (chns == 0) { return; }
    
    int32_t NM = gr->NM;
    for (int32_t km = 0; km < NM; km++)
      { mkgr_mark_t *mk = &(gr->mark.e[km]);
        if (mk->rad > 0.0)
          { double rad = mk->rad * scale;
            bool_t cross = mk->cross;
            bool_t fill = (mk->lwd <= 0.0);
            double xc = mk->ctr.c[0] * scale;
            double yc = mk->ctr.c[1] * scale;
            if (org != NULL) { xc += org->c[0]; yc += org->c[1]; }
            frgb_t color = mk->color;
            double lwd = mk->lwd * scale;
            double ang = mk->ang;
            for (int32_t ic = 0; ic < chns; ic++)
              { /* Get the channel index {c} and the respective mark color {val}: */
                int32_t c = (ch == NULL ? ic : ch[ic]);
                demand((c >= 0) && (c < NC), "invalid channel index");
                float val = color.c[ic % 3];
                /* Get the fill and draw colors: */
                float vfill = (fill ? val : NAN);
                float vdraw = (fill ? NAN : val);
                /* Paint the mark: */
                if (cross)
                  { double ared = fmod(ang, 0.25); /* Angle mod 90 degrees. */
                    bool_t ety = FALSE;
                    bool_t diag = ((ared > 0.0625) && (ared < 0.1875)); 
                    (void)float_image_paint_cross(img, c, xc, yc, rad, ety, lwd/2, diag, vdraw, m);
                  }
                else
                  { bool_t round = TRUE;
                    bool_t diag = FALSE;
                    (void)float_image_paint_dot(img, c, xc, yc, rad, lwd/2, round, diag, vfill, vdraw, m);
                  }
              }
          }
      }
  }

