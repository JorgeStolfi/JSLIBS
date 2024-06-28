/* epswr_vis.h -- window visibility culling procs (in Device coordinates). */
/* Last edited on 2024-06-22 18:37:25 by stolfi */

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

bool_t epswr_vis_segment_is_invisible
  ( epswr_figure_t *epsf,
    double xa, double ya, 
    double xb, double yb
  );

bool_t epswr_vis_curve_is_invisible
  ( epswr_figure_t *epsf,
    double xa, double ya, 
    double xb, double yb, 
    double xc, double yc, 
    double xd, double yd
  );

bool_t epswr_vis_rectangle_is_invisible
  ( epswr_figure_t *epsf,
    double xlo, double xhi, 
    double ylo, double yhi
  );
  
bool_t epswr_vis_centered_rectangle_is_invisible
  ( epswr_figure_t *epsf,
    double xc, double yc, 
    double wd, double ht,
    double ca, double sa
  );
  
bool_t epswr_vis_triangle_is_invisible
  ( epswr_figure_t *epsf,
    double xa, double ya, 
    double xb, double yb, 
    double xc, double yc
  );

bool_t epswr_vis_polygon_is_invisible
  ( epswr_figure_t *epsf,
    double x[], double y[], int32_t npoints
  );

bool_t epswr_vis_circle_is_invisible
  ( epswr_figure_t *epsf,
    double xc, double yc, double rad
  );

bool_t epswr_vis_lune_is_invisible
  ( epswr_figure_t *epsf,
    double xc, double yc, double rad, 
    double tilt
  );

bool_t epswr_vis_slice_is_invisible
  ( epswr_figure_t *epsf,
    double xc, double yc, 
    double rad, 
    double start, double stop
  );
  
#endif
