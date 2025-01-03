#ifndef minn_plot_H
#define minn_plot_H

/* Plots a minimizer's goal function in a 1D or 2D subspace of its domain. */
/* Last edited on 2024-12-05 13:14:30 by stolfi */

#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <float_image.h>

#include <minn.h>
   
void minn_plot_1D_gnuplot
  ( FILE *wr, 
    uint32_t nx, 
    double org[], 
    double u[],
    double rad,
    double step, 
    minn_goal_t *F
  );
  /* Takes a point {org[0..nx-1]} of {\RR^{nx}} and a unit direction
    vector {u[0..nx-1]} in {\RR^{nx}}. Writes to {wr} the values
    {F(nx,org+e*u)} for and odd number of values of {e} spanning {[-rad
    _ +rad]} with the given {step}.
    
    The dimension {nx} must be at least 1. If {org} is {NULL}, assumes
    the origin {(0,..0)}. */

void minn_plot_2D_gnuplot
  ( FILE *wr, 
    uint32_t nx, 
    double org[], 
    double u0[],
    double rad0,
    double u1[],
    double rad1, 
    bool_t box,
    double step,
    minn_goal_t *F
  );
  /* Writes to file {wr} a data for a 3D plot of {F}, over a
    two-dimensional plot domain {\RD} of {\RR^nx} defined by the major
    directions {u0,u1} (unit {nx}-vectors) and the respective radii
    {rad0,rad1}.
    
    The dimension {nx} must be at least 2, and the directions {u0,u1}
    should be orthogonal. The domain {\RD} consisting points of the form
    {org + s0*rad0*u0 + s1*rad1*u1} for various values {s0,s1}. If {org}
    is {NULL}, assumes it is the origin {(0,..0)}.
    
    If if {box} is true, {s0} and {s1} each span the interval {[-1 _ +1]}.
    So {\RD} is a rectangle with sides parallel to {u0} and {u1}
    and half-extents {rad0,rad1} along those directions.
    
    If {box} is false, the pair {(s0,s1)} spans the unit disk of {\RR^2}.
    So {\RD} is the ellipse with major axes parallel to {u0} and {u1}
    and radii {rad0} and {rad1} along those axes.
    
    Either way, the procedure enumerates a grid of points with increment
    {step} on the plane defined by {u0} and {u1}. For each sample point
    {(e0*u0,e1*u1)} of that grid that lies inside the plot domain {\RD},
    it computes the corresponding vector {v = org + e0*u0 + e1*u0} and
    writes to {wr} a line with {e0}, {e1}, and the goal function
    {F(nx,v)}.*/
 
float_image_t *minn_plot_2D_float_image
  ( uint32_t nx, 
    double org[], 
    double u0[],
    double rad0,
    double u1[],
    double rad1, 
    bool_t box,
    double step,
    minn_goal_t *F
  );
  /* Like {minn_plot_2D_gnuplot}, but stores each function value as a
    sample of a {float_image_t}, with {org} at the center, {u0}
    horizontal, and {u1} vertical; then returns that image.
    
    If {box} is false, image samples whose points lie outside the plot 
    domain {\RD} are set to {NAN}. */

#endif
