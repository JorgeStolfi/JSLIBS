/* See cpk_eps.h */
/* Last edited on 2025-01-01 02:44:44 by stolfi */ 

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

#include <cpk_basic.h>

#include <cpk_eps.h>

void cpk_eps_fill_polygon(epswr_figure_t *eps, r2_vec_t *p)
  {
    uint32_t np = p->ne;
    double x[np+2], y[np+2];
    /* For filling, we  copy all vertices from {p} into {x,y},
     inserting the origin (0,0) between every two loops,
     as well as before the first loop and after the last loop.
     This should give the right effect with even/odd filling rule. */

    x[0] = 0; y[0] = 0;
    for (int32_t ip = 0; ip < np; ip++)
      { r2_t *pi = &(p->e[ip]);
        if (r2_is_finite(pi))
          { x[ip+1] = X(*pi); y[ip+1] = Y(*pi); }
        else
          { x[ip+1] = 0; y[ip+1] = 0; }
      }
    x[np+1] = 0; y[np+1] = 0;
    epswr_polygon(eps, TRUE, x, y, (int32_t)np+2, TRUE, FALSE, TRUE);
  }

void cpk_eps_plot_rect_stroke
  ( epswr_figure_t *eps,        /* Where to plot. */
    r2_t *a,             /* Start of stroke. */
    r2_t *b,             /* End of stroke. */
    double r,            /* Half-width of stroke. */
    bool_t fill, bool_t draw /* Plotting options. */
  )
  { 
    /* Zero-area figures contain no pixels: */
    if (r == 0) { return; }
    /* Get the displacement vector from {a} to {b}: */
    r2_t d; r2_sub(b, a, &d);
    /* Zero-area figures contain no pixels: */
    if ((X(d) == 0) && (Y(d) == 0)) { return; }
    /* Get its direction: */
    r2_dir(&d, &d);
    /* Get a vector {v} with length {r} perpendicular to {d}: */
    r2_t v = (r2_t){{-r*Y(d), r*X(d)}};
    /* Build the polygon: */
    double x[5], y[5];
    r2_t p;
    r2_add(a, &v, &p); x[0] = X(p); y[0] = Y(p);
    r2_add(b, &v, &p); x[1] = X(p); y[1] = Y(p);
    r2_sub(b, &v, &p); x[2] = X(p); y[2] = Y(p);
    r2_sub(a, &v, &p); x[3] = X(p); y[3] = Y(p);
    /* Close it: */
    x[4] = x[0]; y[4] = y[0];
    /* Plot the polygon: */
    epswr_polygon(eps, TRUE, x, y, 5, fill, draw, FALSE);
  }

void cpk_eps_plot_circle
  ( epswr_figure_t *eps,        /* Where to plot. */
    r2_t *c,             /* Center of circle. */
    double r,            /* Radius of circle. */
    bool_t fill, bool_t draw /* Plotting options. */
  )
  { /* This one is easy: */
    epswr_circle(eps, X(*c), Y(*c), r, fill, draw);
  }

void cpk_eps_plot_dot
  ( epswr_figure_t *eps,        /* Where to plot. */
    r2_t *c,             /* Center of circle. */
    double r,            /* Radius of circle (mm). */
    bool_t fill, bool_t draw /* Plotting options. */
  )
  { /* This one is easy too: */
    epswr_dot(eps, X(*c), Y(*c), r, fill, draw);
  }

void cpk_eps_draw_polyline(epswr_figure_t *eps, r2_vec_t *p)
  { uint32_t np = p->ne;
    /* We must stroke the segments separately. */
    /* We trust that joints are rounded. */
    double xa = INF, ya = INF;
    for (int32_t ip = 0; ip < np; ip++)
      { r2_t *pb = &(p->e[ip]);
        if (r2_is_finite(pb))
          { double xb = X(*pb);
            double yb = Y(*pb);
            if (xa == INF) 
              { /* Start each polyline with an infinitesimal stroke,
                  so that if it has only one vertex, it will still plot
                  as a dot. (Would a zero-length stroke work too?) */
                xa = xb - 1.0e-6; ya = yb - 1.0e-6;
              }
            epswr_segment(eps, xa, ya, xb, yb);
            xa = xb; ya = yb;
          }
        else
          { xa = INF; ya = INF; }
      }
  }

epswr_figure_t *cpk_eps_new_figure(interval_t B[],  char *fname)
  { 
    double mrg = 2*epswr_pt_per_mm;

    double hPlotSize = 210*epswr_pt_per_mm - 2*mrg;
    double vPlotSize = 297*epswr_pt_per_mm - 2*mrg;
    
    epswr_figure_t *eps = epswr_new_named_figure
      ( NULL, NULL, fname, -1, NULL, 
        hPlotSize, vPlotSize, mrg, mrg, mrg, mrg, TRUE 
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

void cpk_eps_fill_circles
  ( epswr_figure_t *eps,  /* Where to plot. */
    r2_vec_t *c,          /* Centers of circles. */
    double_vec_t *h,      /* Variable radius. */
    double r              /* Fixed radius. */
  )
  { for (int32_t i = 0; i < c->ne; i++)
      { r2_t *ci = &(c->e[i]);
        double hi = ((h == NULL) || (i > h->ne) ? 0 : h->e[i]);
        epswr_circle(eps, X(*ci), Y(*ci), hi+r, TRUE, FALSE);
      }
  }

void cpk_eps_fill_dots
  ( epswr_figure_t *eps,        /* Where to plot. */
    r2_vec_t *c,          /* Centers of circles. */
    double r             /* Dot radius (mm). */
  )
  { for (int32_t i = 0; i < c->ne; i++)
      { r2_t *ci = &(c->e[i]);
        epswr_dot(eps, X(*ci), Y(*ci), r, TRUE, FALSE);
      }
  }

void cpk_eps_show_labels
  ( epswr_figure_t *eps, /* Where to plot. */
    r2_vec_t *c,         /* Reference points. */
    r2_t disp,           /* Extra displacement for all ref pts. */
    double xalign,       /* Horiz rel pos of label to align with ref X. */
    double yalign,       /* Vert rel pos of label to align with ref Y. */
    double ptsize        /* Text point size. */
  )
  { epswr_set_label_font(eps, "Courier", ptsize);
    for (int32_t i = 0; i < c->ne; i++)
      { r2_t *ci = &(c->e[i]);
        r2_t p;
        r2_add(ci, &disp, &p);
        char *txt = jsprintf("%d", i);
        epswr_label(eps, txt, "0", X(p), Y(p), 0.0, FALSE, xalign, yalign, TRUE, FALSE);
        free(txt);
      }
  }

void cpk_eps_draw_circles
  ( epswr_figure_t *eps,        /* Where to plot. */
    r2_vec_t *c,          /* Centers of circles. */
    double_vec_t *h,      /* Variable radius. */
    double r             /* Fixed radius. */
  )
  { for (int32_t i = 0; i < c->ne; i++)
      { r2_t *ci = &(c->e[i]);
        double hi = ((h == NULL) || (i > h->ne) ? 0 : h->e[i]);
        epswr_circle(eps, X(*ci), Y(*ci), hi+r, FALSE, TRUE);
      }
  }
  
void cpk_eps_fill_sausage
  ( epswr_figure_t *eps,        /* Where to plot. */
    r2_vec_t *p,          /* Vertices of sausage. */
    double r             /* Half-width of sausage. */
  )
  { uint32_t np = p->ne;
    r2_t *pa = NULL;
    for (int32_t ip = 0; ip < np; ip++)
      { r2_t *pb = &(p->e[ip]);
        if (r2_is_finite(pb))
          { if (pa != NULL) 
              { cpk_eps_plot_rect_stroke(eps, pa, pb, r, TRUE, FALSE); }
            cpk_eps_plot_circle(eps, pb, r, TRUE, FALSE);
            pa = pb;
          }
        else
          { pa = NULL; }
      }
  }
