/* See {mkgr_mark_grid_draw_eps.h} */
/* Last edited on 2020-11-29 20:32:50 by jstolfi */

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
#include <epswr.h>

#include <mkgr_mark.h>
#include <mkgr_mark_grid.h>
#include <mkgr_mark_grid_draw_eps.h>

void mkgr_mark_grid_draw_eps
  ( epswr_figure_t *eps, 
    mkgr_mark_grid_t *gr, 
    r2_t *org,
    double dpi
  )
  {
    int32_t NM = gr->NM;
    for (int32_t km = 0; km < NM; km++)
      { mkgr_mark_t *mk = &(gr->mark.e[km]);
        double rad = mk->rad;
        if (rad > 0.0)
          { bool_t cross = mk->cross;
            bool_t fill = (mk->lwd == 0.0); /* True if filled only, false if stroked only. */
            bool_t draw = !fill;
            double xc = mk->ctr.c[0];
            double yc = mk->ctr.c[1];
            if (org != NULL) { xc += org->c[0]; yc += org->c[1]; }
            xc = epswr_round_dim_mm(xc, dpi);
            yc = epswr_round_dim_mm(yc, dpi);
            frgb_t color = mk->color;
            double lwd = mk->lwd;
            double ang = mk->ang;
            if (fill)
              { epswr_set_fill_color(eps,  color.c[0], color.c[1], color.c[2]); }
            else
              { epswr_set_pen(eps, color.c[0], color.c[1], color.c[2],  lwd,  0.0, 0.0); }
            if (cross)
              { assert(draw); /* No sense in filling a cross. */
                double ared = fmod(ang, 0.25); /* Angle mod 90 degrees. */
                bool_t diag = ((ared > 0.0625) && (ared < 0.1875)); 
                epswr_cross(eps, xc, yc, rad, diag, draw); 
              }
            else
              { epswr_circle(eps, xc, yc, rad, fill, draw); }
          }
      }
  }

