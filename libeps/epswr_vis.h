/* epswr_vis.h -- window visibility culling procs (in Device coordinates). */
/* Last edited on 2022-10-20 06:51:01 by stolfi */

#ifndef epswr_vis_H
#define epswr_vis_H

#include <epswr.h>
#include <stdint.h>
#include <bool.h>

#include <stdio.h>

/* These procedures return {TRUE} if the figure described by the arguments
  lies certainly outside the current clip path. They return {FALSE}
  when the figure *may* be visible.
  
  All coordinates are in the Device coordinate system; all
  lengths are in pt. */

bool_t epswr_segment_is_invisible
  ( epswr_figure_t *epsf,
    double xa, double ya, 
    double xb, double yb
  );

bool_t epswr_curve_is_invisible
  ( epswr_figure_t *epsf,
    double xa, double ya, 
    double xb, double yb, 
    double xc, double yc, 
    double xd, double yd
  );

bool_t epswr_rectangle_is_invisible
  ( epswr_figure_t *epsf,
    double xlo, double xhi, 
    double ylo, double yhi
  );

bool_t epswr_triangle_is_invisible
  ( epswr_figure_t *epsf,
    double xa, double ya, 
    double xb, double yb, 
    double xc, double yc
  );

bool_t epswr_polygon_is_invisible
  ( epswr_figure_t *epsf,
    double x[], double y[], int32_t npoints
  );

bool_t epswr_circle_is_invisible
  ( epswr_figure_t *epsf,
    double xc, double yc, double rad
  );

bool_t epswr_lune_is_invisible
  ( epswr_figure_t *epsf,
    double xc, double yc, double rad, 
    double tilt
  );

bool_t epswr_slice_is_invisible
  ( epswr_figure_t *epsf,
    double xc, double yc, 
    double rad, 
    double start, double stop
  );
  
#endif
