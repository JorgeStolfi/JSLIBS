/* See fltgraph.h */
/* Last edited on 2021-06-26 18:39:56 by jstolfi */

#include <fltgraph.h>
#include <affirm.h>
#include <bool.h>
#include <flt.h>
#include <epswr.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void fltgraph_plot
  ( epswr_figure_t *fig,
    Float f (Float x),
    Interval xd,
    Interval yd,
    int m
  )
  { int xi;
    Float x0, y0, x1, y1;

    ROUND_NEAR;
    epswr_comment(fig, "Plot of actual graph");

    x1 = xd.lo;
    y1 = f(x1);

    for (xi=0; xi<m; xi++)
      { x0 = x1;
        y0 = y1;

        x1 = (Float)(xd.lo + ((xd.hi - xd.lo)*(double)(xi+1))/(double)m);
        y1 = f(x1);

        if 
          ( ((fabs(y1) < Infinity) && (fabs(y0) < Infinity)) &&
            ((y0 >= yd.lo) || (y1 >= yd.lo)) &&
            ((y0 <= yd.hi) || (y1 <= yd.hi))
          )
          { epswr_segment(fig, x0, y0, x1, y1); }
      }
  }

void fltgraph_draw_axes
  ( epswr_figure_t *fig,
    Interval xd,
    Interval yd
  )
  { if ((xd.lo < Zero) && (xd.hi > Zero))
      { epswr_coord_line(fig, epswr_axis_HOR, 0.0); }
    if ((yd.lo < Zero) && (yd.hi > Zero))
      { epswr_coord_line(fig, epswr_axis_VER, 0.0); }
  }

#define MAXLABLEN (300)

void fltgraph_draw_tics
  ( epswr_figure_t *fig,
    epswr_axis_t axis,
    Float lo, Float hi,
    int n,
    double ticsz,
    char *labfmt,
    double labalign,
    Interval clip
  )
  { 
    int i;
    char buf[MAXLABLEN];
    double eps = 0.1*(hi - lo)/((double)n); /* Origin-avoidance radius. */
    for (i = 0; i <= n; i++)
      { double r = ((double)i)/((double)n);
        double c = (1-r)*lo + r*hi;
        if ((c > clip.lo) && (c <= clip.hi) && (fabs(c) > eps))
          { if (labfmt != NULL) { snprintf(buf, MAXLABLEN, labfmt, c); }
            if (axis == 0)
              { epswr_tic(fig, axis, c, 0.0, ticsz, 0.5);
                if (labfmt != NULL) 
                  { epswr_label(fig, buf, "0", c,0.0, 0.0, TRUE, 0.5,labalign, TRUE,FALSE); }
              }
            else if (axis == 1) 
              { epswr_tic(fig, axis, 0.0, c, ticsz, 0.5);
                if (labfmt != NULL) 
                  { epswr_label(fig, buf, "0", 0.0,c, 0.0, TRUE, labalign,0.5, TRUE,FALSE); }
              }
            else
              { affirm(FALSE, "bad axis"); }
          }
      }
    epswr_flush(fig);
  }    
