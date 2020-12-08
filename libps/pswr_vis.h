/* pswr_vis.h -- window visibility culling procs (in device coordinates). */
/* Last edited on 2009-01-06 03:42:02 by stolfi */

#ifndef pswr_vis_H
#define pswr_vis_H

#include <pswr.h>
#include <bool.h>

#include <stdio.h>

/* These procedures return {TRUE} if the figure described by the arguments
  lies certainly outside the current clip path. They return {FALSE}
  when the figure *may* be visible.
  
  All coordinates are in the device coordinate system; all
  lengths are in pt. */

bool_t pswr_segment_is_invisible
  ( PSStream *ps,
    double xa, double ya, 
    double xb, double yb
  );

bool_t pswr_curve_is_invisible
  ( PSStream *ps,
    double xa, double ya, 
    double xb, double yb, 
    double xc, double yc, 
    double xd, double yd
  );

bool_t pswr_rectangle_is_invisible
  ( PSStream *ps,
    double xlo, double xhi, 
    double ylo, double yhi
  );

bool_t pswr_triangle_is_invisible
  ( PSStream *ps,
    double xa, double ya, 
    double xb, double yb, 
    double xc, double yc
  );

bool_t pswr_polygon_is_invisible
  ( PSStream *ps,
    double x[], double y[], int npoints
  );

bool_t pswr_circle_is_invisible
  ( PSStream *ps,
    double xc, double yc, double rad
  );

bool_t pswr_lune_is_invisible
  ( PSStream *ps,
    double xc, double yc, double rad, 
    double tilt
  );

bool_t pswr_slice_is_invisible
  ( PSStream *ps,
    double xc, double yc, 
    double rad, 
    double start, double stop
  );
  
#endif
