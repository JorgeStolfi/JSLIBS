#ifndef mkgr_mark_H
#define mkgr_mark_H

/* mkgr_mark.h - parameters that describe a mark in a grid. */
/* Last edited on 2020-11-29 17:36:34 by jstolfi */

#define _GNU_SOURCE
#include <bool.h>
#include <vec.h>
#include <r2.h>
#include <frgb.h>

typedef struct mkgr_mark_t
  { r2_t ctr;     /* Center coordinates of each mark (mm). */
    bool_t cross; /* Type of each mark: TRUE if cross, FALSE if circle. */
    double rad;   /* Radius of mark (mm). */
    double ang;   /* Rotation angle of each mark (as fraction of full turn). */
    frgb_t color; /* Color of mark. */
    double lwd;   /* Line width for stroking (mm). */
  } mkgr_mark_t;
  /* Data for a mark in a grid, centered at {ctr}.
  
    If {cross} is true, the mark is a cross.  Then {lwd} must
    be positive.  The radius {rad} is the extent of the arms,
    not including the pen radius;  so the actual radius is {rad+lwd/2}.
    
    If {cross} is false, the mark is a disk or circle
    with radius {rad}.  The angle {ang} must be zero.
    In this case, if {lwd} is zero, the mark is a disk filled with the given {color},
    with no outline.  If {lwd} is positive, the mark is an unfilled
    circle stroked in the given {color} with a pen of width {lwd},
    and its actual radius will be {rad+lwd/2}. */
  
mkgr_mark_t mkgr_make_dot(r2_t ctr, double rad, frgb_t color);
  /* Creates a mark that is a disk with center {ctr}, radius {rad},
    filled with the given {color}, without stroked border. */

mkgr_mark_t mkgr_make_circle(r2_t ctr, double rad, double lwd, frgb_t color);
  /* Creates a mark that is an unfilled circle with center {ctr}, radius {rad},
    stroked with the given {color} and line width {lwd}. */

mkgr_mark_t mkgr_make_cross(r2_t ctr, double rad, double ang, double lwd, frgb_t color);
  /* Creates a mark that is a cross with center {ctr}, nominal radius {rad},
    rotated by {ang} turns, stroked with the given {color} and line width {lwd}. */

vec_typedef(mkgr_mark_vec_t, mkgr_mark_vec, mkgr_mark_t);
  /* An extensible vector of marks. */

#endif
