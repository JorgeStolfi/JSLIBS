#ifndef float_image_interpolate_H
#define float_image_interpolate_H

/* Bilinear (C0) and Bicubic (C1) interpolation of floating-point images. */
/* Last edited on 2012-01-07 00:02:33 by stolfilocal */ 

#include <bool.h>
#include <indexing.h>
#include <float_image.h>

/* 
  PIXEL INDICES AND COORDINATES
  
  The procedures in this module assume that the pixel on column {ix} and 
  row {iy} is a unit square centered at point {(ix+0.5,iy+0.5)}.  
  
  PIXELS OUTSIDE THE IMAGE DOMAIN
  
  Pixel indices outside the valid ranges {0..NX-1}, {0..NY-1} are
  reduced according to {ix_reduce(i, N, red)}. If a reduced pixel
  index is {-1}, the functions assume that its sample values are all 0.
  
  INTERPOLATION METHODS
  
  If the {order} parameter is -1, the interpolation is piecewise
  constant and discontinuous: the result is the value of the pixel
  that contains the given point {(x,y)}.
  
  If {order} is 0, the interpolation is piecewise biaffine (degree 1)
  and C0. The result is a weighted average of the 2×2 block of pixels whose
  centers surround the point {(x,y)}. It reconstructs bilinear
  functions.
  
  If {order} is 1, the interpolation is bicubic and C1. The result is
  a weighted average of the 4×4 block of pixels such that the point
  {(x,y)} is surrounded by the centers of the innermost 2×2 block of
  pixels. It reconstructs biquadratic functions. */

double float_image_interpolate_sample
  ( float_image_t *A, 
    int c, 
    double x, 
    double y, 
    int order, 
    ix_reduction_t red
  );
  /* Returns the value of channel {c} of image {A} at the point
    {(x,y)}, computed by interpolation from nearby samples. */

void float_image_interpolate_pixel
  ( float_image_t *A, 
    double x, 
    double y, 
    int order, 
    ix_reduction_t red, 
    double z[]
  );
  /* For {c=0..NC-1}, stores into {z[c]} the value of chanel {c} of
    image {A} at the point {(x,y)}, computed as in
    {float_image_interpolate_sample}; where {NC} is the number of
    channels of {A}. */

/* 
  INTERPOLATION ON A REGULAR GRID
  
  The procedures in this section refer to a grid of sampling points,
  defined by center coordinates {rx,ry}, non-negative half-extents
  {hx,hy}, and coordinate increments {dx,dy}. Namely, the sampling
  points are {p[-hx..+hx,-hy..+hy]} where {p[jx,jy]} is the point
  {(rx+jx*dx, ry+jy*dy)}.  Thus the grid has {nx=2*hx+1}
  columns and {ny=2*hy+1} rows.  */

void float_image_interpolate_grid_samples
  ( float_image_t *A, 
    int c, 
    double rx, int hx, double dx,
    double ry, int hy, double dy,
    int order, 
    ix_reduction_t red,
    double z[]
  );
  /* 
    Interpolates channel {c} of image {A} at the sampling point grid
    defined by parameters {rx,hx,dx,ry,hy,dy}.  

    The interpolated sample values are stored by rows as a single vector {z[0..nx*ny-1]}, 
    where {nx = 2*hx+1} and {ny = 2*hy+1}; so that the interpolated value at {p[jx,jy]} is 
    {z[jxy]} where {jxy = (jx+hx)+(jy+hy)*nx}. */

void float_image_interpolate_grid_pixels
  ( float_image_t *A, 
    double rx, int hx, double dx,
    double ry, int hy, double dy,
    int order, 
    ix_reduction_t red, 
    double z[]
  );
  /* 
    Applies {float_image_interpolate_grid_samples} to each channel 
    of {A}, with parameters {x,hx,dx,y,hy,dy,order,red}. 
    The interpolated channel {c} value at point {p[x,y]}
    is stored in {z[NC*jxy + c]} where {jxy = (jx+hx)+(jy+hy)*nx}. */


#endif
