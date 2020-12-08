/* See pswr_vis.h */
/* Last edited on 2008-05-25 01:13:31 by stolfi */

#define _GNU_SOURCE
#include <pswr.h>
#include <pswr_vis.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>

#define Sqrt3 1.73205080756887729352

#define DEFWINDOW \
    double xMin = ps->hMin - 3.0; \
    double xMax = ps->hMax + 3.0; \
    double yMin = ps->vMin - 3.0; \
    double yMax = ps->vMax + 3.0

#define INITBBOX \
    double xlo = +INFINITY; \
    double xhi = -INFINITY; \
    double ylo = +INFINITY; \
    double yhi = -INFINITY
    
#define XBOX(x) { if ((x) < xlo) { xlo = (x); } if ((x) > xhi) { xhi = (x); } }
#define YBOX(y) { if ((y) < ylo) { ylo = (y); } if ((y) > yhi) { yhi = (y); } }

#define CHECKBOX \
    return ((xlo > xMax) || (xhi < xMin) || (ylo > yMax) || (yhi < yMin))

bool_t pswr_segment_is_invisible
  ( PSStream *ps,
    double xa, double ya, 
    double xb, double yb
  )
  { DEFWINDOW;
    INITBBOX;
    XBOX(xa); YBOX(ya);
    XBOX(xb); YBOX(yb);
    CHECKBOX;
  }

bool_t pswr_curve_is_invisible
  ( PSStream *ps,
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

bool_t pswr_rectangle_is_invisible
  ( PSStream *ps,
    double xlo, double xhi, 
    double ylo, double yhi
  )
  { DEFWINDOW;
    CHECKBOX;
  }

bool_t pswr_triangle_is_invisible
  ( PSStream *ps,
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

bool_t pswr_polygon_is_invisible
  ( PSStream *ps,
    double x[], double y[], int npoints
  )
  { DEFWINDOW;
    INITBBOX;
    int i;
    for (i = 0; i < npoints; i++) { XBOX(x[i]); YBOX(y[i]); }
    CHECKBOX;
  }

bool_t pswr_circle_is_invisible
  ( PSStream *ps,
    double xc, double yc, double rad
  )
  { DEFWINDOW;
    INITBBOX;
    XBOX(xc-rad); YBOX(yc-rad);
    XBOX(xc+rad); YBOX(yc+rad);
    CHECKBOX;
  }


bool_t pswr_lune_is_invisible
  ( PSStream *ps,
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

bool_t pswr_slice_is_invisible
  ( PSStream *ps,
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
