/* See fgraph.h */

#include "fgraph.h"
#include "foifloat.h"
#include "foimisc.h"
#include "iomisc.h"
#include "interval.h"
#include "plt0.h"
#include <math.h>
#include <stdio.h>

void fgraph_plot(
    FILE *psfile,
    Float f (Float x),
    Interval xd,
    Interval yd,
    int m
  )
  {
    int xi;
    Float x0, y0, x1, y1;

    ROUND_NEAR;
    plt0_begin_section(psfile, "Plot of actual graph");

    plt0_set_pen(psfile, 0.0, 0.30, 0.0, 0.0);

    x1 = xd.lo;
    y1 = f(x1);

    for (xi=0; xi<m; xi++)
      {
        x0 = x1;
        y0 = y1;

        x1 = xd.lo + ((xd.hi - xd.lo)*(xi+1))/m;
        y1 = f(x1);

        if ( ((abs(y1) < Infinity) && (abs(y0) < Infinity))
          && ((y0 <= yd.hi) || (y1 <= yd.hi))
          && ((y0 >= yd.lo) || (y1 >= yd.lo))
          )
          { plt0_draw_segment(psfile, x0, y0, x1, y1); }

      }

    plt0_end_section(psfile);
  }

void fgraph_draw_axes(
    FILE *psfile,
    Interval xd,
    Interval yd
  )
  {
    if ((xd.lo < Zero) && (xd.hi > Zero))
      { plt0_draw_coord_line(psfile, 'x', 0.0); }
    if ((yd.lo < Zero) && (yd.hi > Zero))
      { plt0_draw_coord_line(psfile, 'y', 0.0); }
  }

