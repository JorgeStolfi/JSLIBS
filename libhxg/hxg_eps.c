/* See hxg_eps.h */
/* Last edited on 2025-01-01 02:35:11 by stolfi */ 

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <r2.h>
#include <i2.h>
#include <interval.h>
#include <jsprintf.h>
#include <epswr.h>
#include <frgb.h>

#include <cpk_basic.h>
#include <hxg_canvas.h>
#include <hxg_paint.h>
#include <hxg_eps.h>

epswr_figure_t *hxg_eps_new_figure(interval_t B[],  char *fname)
  { 
    double mrg = 2*epswr_pt_per_mm;

    double hPlotSize = 210*epswr_pt_per_mm - 2*mrg;
    double vPlotSize = 297*epswr_pt_per_mm - 2*mrg;
    
    bool_t verbose = FALSE;
    epswr_figure_t *eps = epswr_new_named_figure
      ( NULL, NULL, fname, -1, NULL, 
        hPlotSize, vPlotSize, mrg, mrg, mrg, mrg, verbose 
      );
    
    /* Compute the client plotting extents and ranges with some slop: */
    double zxraw = HI(B[0]) - LO(B[0]);
    double zyraw = HI(B[1]) - LO(B[1]);
    double mx = 0.05 * zxraw, my = 0.05 * zyraw; /* Margins */
    double xmin = LO(B[0]) - mx, xmax = HI(B[0]) + mx;
    double ymin = LO(B[1]) - my, ymax = HI(B[1]) + my;
    epswr_set_client_window(eps, xmin, xmax, ymin, ymax);

    return eps;
  }

void hxg_eps_plot_canvas
  ( epswr_figure_t *eps,  /* Where to plot. */
    hxg_canvas_t *cvs,    /* The canvas to plot. */
    double pixv[],        /* Pixel values, stored in the standard way. */
    uint32_t maxv,        /* Max pixel value in tables. */
    double *r,            /* Radius table. */
    frgb_t *fill_color,  /* Fill color table. */
    frgb_t *draw_color   /* Draw color table. */
  )
  {
    bool_t debug = FALSE;
    /* Get hold of params: */
    uint32_t nx = cvs->size[0]; /* Num of points per canvas row. */
    uint32_t ny = cvs->size[1]; /* Num of rows in canvas. */
    /* Current pixel value, initial state not relevant: */
    double cv = 0;
    /* Current dot radius: */
    double cr = INF;
    /* Current Postscript state, starts undefined: */
    frgb_t cf = (frgb_t){{NAN,NAN,NAN}}; /* Current fill color. */
    frgb_t cd = (frgb_t){{NAN,NAN,NAN}}; /* Current draw color. */
    /* Plot pixels: */
    for (uint32_t iy = 0; iy < ny; iy++)
      { uint32_t k = iy*nx;
        for (uint32_t ix = 0; ix < nx; ix++,k++)
          { double val = pixv[k];
            if ((cr == INF) || (val != cv))
              { /* Round {val} to integer: */
                int32_t iv = (fabs(val) > (double)maxv ? (int32_t)maxv : (int32_t)floor(fabs(val) + 0.5));
                if (val < 0) { iv = -iv; }
                if (debug) { fprintf(stderr, "  val = %24.16e  iv = %+d\n", val, iv); }
                /* Get radius: */
                cr = r[iv];
                if (cr > 0) 
                  { cf = fill_color[iv];
                    cd = draw_color[iv];
                    if (val < 0)
                      { /* Complement {cf}: */
                        if ((! isfinite(cf.c[0])) || (cf.c[0] < 0))
                          { cf = (frgb_t){{0,0,0}}; }
                        else
                          { cf = (frgb_t){{ 1-cf.c[0], 1-cf.c[1], 1-cf.c[2] }}; }
                      }
                    epswr_set_pen(eps, cd.c[0], cd.c[1], cd.c[2], 0.1, 0.0, 0.0);
                    epswr_set_fill_color(eps, cf.c[0], cf.c[1], cf.c[2]);
                  }
                cv = val;
              }
            if (cr > 0) 
              { r2_t c = hxg_canvas_pixel_pos(cvs, (i2_t){{ (int32_t)ix, (int32_t)iy }});
                epswr_dot(eps, X(c), Y(c), cr, TRUE, TRUE);
              }
          }
      }
  }
