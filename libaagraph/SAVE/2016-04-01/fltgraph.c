/* See fltgraph.h */
/* Last edited on 2007-12-26 19:17:40 by stolfi */

#include <fltgraph.h>
#include <affirm.h>
#include <bool.h>
#include <flt.h>
#include <pswr.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void fltgraph_plot
  ( PSStream *ps,
    Float f (Float x),
    Interval xd,
    Interval yd,
    int m
  )
  { int xi;
    Float x0, y0, x1, y1;

    ROUND_NEAR;
    pswr_comment(ps, "Plot of actual graph");

    x1 = xd.lo;
    y1 = f(x1);

    for (xi=0; xi<m; xi++)
      { x0 = x1;
        y0 = y1;

        x1 = xd.lo + ((xd.hi - xd.lo)*(xi+1))/m;
        y1 = f(x1);

        if 
          ( ((abs(y1) < Infinity) && (abs(y0) < Infinity)) &&
            ((y0 >= yd.lo) || (y1 >= yd.lo)) &&
            ((y0 <= yd.hi) || (y1 <= yd.hi))
          )
          { pswr_segment(ps, x0, y0, x1, y1); }
      }
  }

void fltgraph_draw_axes
  ( PSStream *ps,
    Interval xd,
    Interval yd
  )
  { if ((xd.lo < Zero) && (xd.hi > Zero))
      { pswr_coord_line(ps, HOR, 0.0); }
    if ((yd.lo < Zero) && (yd.hi > Zero))
      { pswr_coord_line(ps, VER, 0.0); }
  }

#define MAXLABLEN (300)

void fltgraph_draw_tics
  ( PSStream *ps,
    pswr_axis_t axis,
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
              { pswr_tic(ps, axis, c, 0.0, ticsz, 0.5);
                if (labfmt != NULL) 
                  { pswr_label(ps, buf, c, 0.0, 0.5, labalign); }
              }
            else if (axis == 1) 
              { pswr_tic(ps, axis, 0.0, c, ticsz, 0.5);
                if (labfmt != NULL) 
                  { pswr_label(ps, buf, 0.0, c, labalign, 0.5); }
              }
            else
              { affirm(FALSE, "bad axis"); }
          }
      }
    pswr_flush(ps);
  }    
