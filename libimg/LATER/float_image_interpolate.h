#ifndef float_image_interpolate_H
#define float_image_interpolate_H

/* Bilinear (C0) and Bicubic (C1) interpolation of floating-point images. */
/* Last edited on 2024-12-04 23:38:20 by stolfi */ 

#include <bool.h>
#include <ix_reduce.h>
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
  
  The interpolation method is defined by the {order} and {convex}
  parameters. In every case the output is a piecewise polinomial
  function of the {x} and {y} parameters, with pieces of degree
  {order+1}. The function depends on a window of {W×W} samples whose
  center is closest to the point {(x,y)}, where {W} depends on {order}
  and {convex}.
  
  The {order} parameter is the continuity order of the interpolated
  function, and defines (indirectly) how many samples of the input image are used
  to compute the interpolated value for a given point.  Possible values
  are {-1} (result is piecewise constant, with step-like discontinuities),
  {0} (bilinear, continuous with derivative discontinities),
  {1} (continuous function and first derivatives) and {2} (continuous function.
  first and second derivatives).
  
  The {convex} flag is relevant only when {order>=1},and selects
  between two interpolation stategies:

    If the {convex} flag is true, the result is always a convex
    combination of the values of the image samples in the window, and
    therefore always lies in the interval spanned by those sample
    values. The window size {W} is {order+2}. Indeed, it is the
    tensor-type B-spline function of degree {order+1} with those pixels
    as control values. However, when {order} is 1 or more, the result
    will not be a true interpolation: when evaluated at pixel centers it
    will return a slightly smoothed version of the original image. The
    derivatives of order {order+1} may be discontinuous at integer
    values of {x} or {y}.

    If the {convex} flag is false, the function is always a true
    interpolation: when evaluated at a pixel center, it returns the
    original value of that pixel. Otherwise the window size {W} is {1} if
    {order=-1}, or {2*(order+1)} when{order>=0}, The result may over- or
    undershoot the range spanned by the samples in the window, by at most
    {???} times the difference between the largest and smallest sample in
    the window. The derivatives of order {order+1} may be discontinous at
    integer and (if {order >= 1}) half-integer values of {x} or {y}.

  Said otherwise:
  
  If {order} is -1, the interpolation is piecewise constant and
  discontinuous: the result is the value of the image pixel that
  contains the given point {(x,y)}. The {convex} flag is irrelevant,
  since this this interpolation is both exact and convex.
  
  If {order} is 0, the interpolation is piecewise biaffine (i.e., of
  degree 1 in each coordinate) and C0. It reconstructs bilinear
  functions. The {convex} flag is irrelevant, since this this
  interpolation is both exact and convex.
  
  If {order} is 1, then the interpolation is piecewise biquadratic and
  C1. It is exact if {convex} is false, convex if {convex} is true.
  
  If {order} is 2, then the interpolation is piecewise bicubic and
  C2. It is exact if {convex} is false, convex if {convex} is true. */

double float_image_interpolate_sample
  ( float_image_t *A, 
    int32_t c, 
    double x, 
    double y, 
    int32_t order, 
    bool_t convex,
    ix_reduce_mode_t red
  );
  /* Returns the value of channel {c} of image {A} at the point
    {(x,y)}, computed by interpolation from nearby samples. */

void float_image_interpolate_pixel
  ( float_image_t *A, 
    double x, 
    double y, 
    int32_t order, 
    bool_t convex,
    ix_reduce_mode_t red, 
    double z[]
  );
  /* For {c=0..NC-1}, stores into {z[c]} the value of chanel {c} of
    image {A} at the point {(x,y)}, computed as in
    {float_image_interpolate_sample}; where {NC} is the number of
    channels of {A}. */

/* 
  MULTIPLE INTERPOLATIONS ON A REGULAR GRID
  
  The procedures in this section refer to a grid of sampling points,
  defined by center coordinates {rx,ry}, non-negative half-extents
  {hx,hy}, and coordinate increments {dx,dy}. Namely, the sampling
  points are {p[-hx..+hx,-hy..+hy]} where {p[jx,jy]} is the point
  {(rx+jx*dx, ry+jy*dy)}.  Thus the grid has {nx=2*hx+1}
  columns and {ny=2*hy+1} rows.  */

void float_image_interpolate_grid_samples
  ( float_image_t *A, 
    int32_t c, 
    double rx, int32_t hx, double dx,
    double ry, int32_t hy, double dy,
    int32_t order, 
    bool_t convex,
    ix_reduce_mode_t red,
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
    double rx, int32_t hx, double dx,
    double ry, int32_t hy, double dy,
    int32_t order, 
    bool_t convex,
    ix_reduce_mode_t red, 
    double z[]
  );
  /* 
    Applies {float_image_interpolate_grid_samples} to each channel 
    of {A}, with parameters {x,hx,dx,y,hy,dy,order,red}. 
    The interpolated channel {c} value at point {p[x,y]}
    is stored in {z[NC*jxy + c]} where {jxy = (jx+hx)+(jy+hy)*nx}. */

#endif
