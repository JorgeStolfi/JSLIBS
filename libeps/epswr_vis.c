/* See epswr_vis.h */
/* Last edited on 2024-12-05 10:15:06 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>

#include <epswr.h>
#include <epswr_def.h>
#include <epswr_vis.h>

#define Sqrt3 1.73205080756887729352

#define DEFWINDOW \
    double xMin = epsf->hMin - 3.0; \
    double xMax = epsf->hMax + 3.0; \
    double yMin = epsf->vMin - 3.0; \
    double yMax = epsf->vMax + 3.0

#define INITBBOX \
    double xlo = +INFINITY; \
    double xhi = -INFINITY; \
    double ylo = +INFINITY; \
    double yhi = -INFINITY
    
#define XBOX(x) { if ((x) < xlo) { xlo = (x); } if ((x) > xhi) { xhi = (x); } }
#define YBOX(y) { if ((y) < ylo) { ylo = (y); } if ((y) > yhi) { yhi = (y); } }

#define CHECKBOX \
    return ((xlo > xMax) || (xhi < xMin) || (ylo > yMax) || (yhi < yMin))

bool_t epswr_vis_segment_is_invisible
  ( epswr_figure_t *epsf,
    double xa, double ya, 
    double xb, double yb
  )
  { DEFWINDOW;
    INITBBOX;
    XBOX(xa); YBOX(ya);
    XBOX(xb); YBOX(yb);
    CHECKBOX;
  }

bool_t epswr_vis_curve_is_invisible
  ( epswr_figure_t *epsf,
    double xa, double ya, 
    double xb, double yb, 
    double xc, double yc, 
    double xd, double yd
  )
  { DEFWINDOW;
    INITBBOX;
    XBOX(xa); YBOX(ya);
    XBOX(xb); YBOX(yb);
    XBOX(xc); YBOX(yc);
    XBOX(xd); YBOX(yd);
    CHECKBOX;
  }

bool_t epswr_vis_rectangle_is_invisible
  ( epswr_figure_t *epsf,
    double xlo, double xhi, 
    double ylo, double yhi
  )
  { DEFWINDOW;
    CHECKBOX;
  }

bool_t epswr_vis_centered_rectangle_is_invisible
  ( epswr_figure_t *epsf,
    double xc, double yc, 
    double wd, double ht,
    double ca, double sa
  )
  { /* Get bounding box coords {psxlo,psxhi,psylo,psyhi}: */
    double xr = (fabs(ca)*wd + fabs(sa)*ht)/2;
    double yr = (fabs(sa)*wd + fabs(ca)*ht)/2;
    double xlo = xc - xr;
    double xhi = xc + xr;
    double ylo = yc - yr;
    double yhi = yc + yr;
    return epswr_vis_rectangle_is_invisible(epsf, xlo, xhi, ylo, yhi);
  }

bool_t epswr_vis_triangle_is_invisible
  ( epswr_figure_t *epsf,
    double xa, double ya, 
    double xb, double yb, 
    double xc, double yc
  )
  { DEFWINDOW;
    INITBBOX;
    XBOX(xa); YBOX(ya);
    XBOX(xb); YBOX(yb);
    XBOX(xc); YBOX(yc);
    CHECKBOX;
  }

bool_t epswr_vis_polygon_is_invisible
  ( epswr_figure_t *epsf,
    double x[], double y[], int32_t npoints
  )
  { DEFWINDOW;
    INITBBOX;
    int32_t i;
    for (i = 0; i < npoints; i++) { XBOX(x[i]); YBOX(y[i]); }
    CHECKBOX;
  }

bool_t epswr_vis_circle_is_invisible
  ( epswr_figure_t *epsf,
    double xc, double yc, double rad
  )
  { DEFWINDOW;
    INITBBOX;
    XBOX(xc-rad); YBOX(yc-rad);
    XBOX(xc+rad); YBOX(yc+rad);
    CHECKBOX;
  }


bool_t epswr_vis_lune_is_invisible
  ( epswr_figure_t *epsf,
    double xc, double yc, double rad, 
    double tilt
  )
  { DEFWINDOW;
    INITBBOX;
    /* The tip of the lune lies "rad*sqrt(3)" away from center: */ 
    double height = rad*Sqrt3;
    XBOX(xc-height); YBOX(yc-height);
    XBOX(xc+height); YBOX(yc+height);
    CHECKBOX;
  }

bool_t epswr_vis_slice_is_invisible
  ( epswr_figure_t *epsf,
    double xc, double yc, 
    double rad, 
    double start, double stop
  )
  { DEFWINDOW;
    INITBBOX;
    /* Use whole circle, igoring "start,stop": */
    XBOX(xc-rad); YBOX(yc-rad);
    XBOX(xc+rad); YBOX(yc+rad);
    CHECKBOX;
  }
