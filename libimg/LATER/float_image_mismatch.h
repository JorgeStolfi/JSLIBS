#ifndef float_image_mismatch_H
#define float_image_mismatch_H

/* Tools for local or global comparison of two or more images. */
/* Last edited on 2024-12-04 23:21:37 by stolfi */ 

#include <r2.h>
#include <i2.h>

/* 
  IMAGE MODEL

  The mismatch functions in this interface assume that each image is a
  noisy real-valued function of two real variables, identified by an /image
  index/ {i} in some given range {0..ni-1}.
  
  More precisely, for each image index {i} in {0..ni-1} and each point
  {p} of the plane, there is a real /image mean value/ {z[i](p)} and a
  corresponding /image variance/ {v[i](p)}. The true value of the
  image at that point is assumed to be {z[i](p) + e[i](p)}, where the
  values {e[i](p)} are unknown variables, independently drawn from a
  Gaussian distribution with zero mean and variance {v[i](p)}.

  If the point {p} falls inside the image domain, the variance
  {v[i](p)} should normally be the variance of the image capture noise
  plus the variance of the quantization error. If the quantization
  error is uniformly distributed in the interval {[-q/2 _ +q/2]}, for
  example, its variance alone should be about {q^2/12}.

  If the sample point falls outside the image domain, the quantities
  {z[i](p)} and {v[i](p)} should be the expected value and variance of
  a completely unknown pixel. If image values are uniformly
  distributed in {[0_1]}, then {z[i](p)} should be {1/2} and {v[i](p)}
  should be {1/12}.

  SAMPLING GRID

  The mismatch functions below are computed by evaluating each image
  on a /grid of sampling points/, defined by a given 
  /reference point/ {r} and two naturals {hx,hy}. 
  
  The grid consists of {nx*ny} points {p[0..nx*ny-1]} of the plane,
  where {nx = 2*hx+1} and {ny = 2*hy+1}. Each sampling point {p[k]} is
  the point {r+(jx,jy)}, where {jx} ranges in {-hx..+hx}, {jy} ranges
  over {-hy..+hy}, and {k = (jx+hx) + nx*(jy+hy)}.
  
  Note that the sampling grid is centered at the given reference point
  {r}. Note also that the sampling points are spaced one unit apart in
  both axes; but since the coordinates of {r} may be arbitrary
  fractional numbers, the sampling points {p[k]} may be fractional
  too. */

typedef void float_image_eval_t(int32_t i, r2_t *r, int32_t hx, int32_t hy, double z[], double v[]); 
  /* Type of a procedure that evaluates image number {i} at the 
    points of the sampling grid defined by the center {r} and
    the integers {hx,hy}. 
    
    The procedure shoud return the image values in {z[0..nx*ny-1]}, and
    the corresponding image variances in {v[0..nx*ny-1]}, row-by-row. */
    
/*  
  MISMATCH FUNCTION
  
  For this interface, an /image mismatch function/ is a function {Q(r)}
  that measures the overall difference betwen {ni} images over the 
  sampling grids defined by the reference points {r[0..ni-1]}
  and some grid half-sizes {hx,hy}.  (The values {ni,hi,xi} are
  also parameters of {Q}, but we omit them for simplicity.) 
  
  The result of an image mismatch function is always non-negative. It
  is zero when the images are identical over all corresponding points
  of the sampling grids. A mismatch function is zero when {ni = 1},
  and is undefined (NaN) when {ni = 0}.
  
  Each of the mismatch functions defined here is /smooth/, meaning
  that the value of {Q(r+t)}, where {t[0..ni-1]} are reasonable
  perturbations, is well approximated by a quadratic function of
  {t[0..ni-1]}. Here /reasonable perturbation/ means that the absolute
  value of each coordinate of each {t[i]} is at most {1/2}, but still
  large enough to make roudoff errors negligible; and /well
  approximated/ means that the residual of the approximation is much
  smaller in absolute value than the terms of order 1 and 2.
  
  !!! What about the triangle inequality? !!!
*/

!!! STOPPED HERE 2009-07-10 01:07:25

  /*
  
  
  WINDOW WEIGHT MASK
  
  An image mismatch function accepts also a /window weight/
  for each point of the sampling grid. The window weights
  are defined by   the grid half-sizes {hx,hy} and a pair of one-dimensional weight tables
  {wx[0..nx-1]} and {wy[0..ny-1]}.  Specifically, 
  the window weight for a sampling point {p[i][k] = r[i]+(jx,jy)} is 
  {W(
  
  {w(p)} for any given point {p}. These weights are used when
  computing averages over several points {p}.
  
  SAMPLING FUNCTION
  
  The images to be compared are defined implicitly by a
  client-specified /sampling function/ {eval} of type {fimm_eval_t}. */

double float_image_mismatch_var
  ( int32_t ni,                   /* Number of images being compared. */
    float_image_eval_t *eval, /* Single-channel image evaluator. */
    int32_t hx,                  /* Half-width of comparison window. */
    double wx[],              /* Horizontal weight table. */
    int32_t hy,                  /* Half-height of comparison window. */
    double wy[],              /* Vertical weight table. */
    r2_t r[]                  /* Sampling grid center for each image. */
  );
  /* Returns a mismatch measure for images {0..ni-1}, based on the
    mean squared differences of corresponding sample values
    This is essentially the weighted mean-square-difference mismatch,
    generalized to {ni} images.
    
    For each image {i}, the procedure generates an array of {nx*ny}
    sampling points {p[i,0..nx*ny-1]} centered on the fractional
    point {r[i]} and spaced one pixel apart; where {nx = 2*hx+1} and
    {ny = 2*hy+1}.
    
    For each {k} in {0..nx*ny-1}, the procedure computes a
    /quadratic point mismatch/ {V2[k]}. The result of the function is
    the weighted average of those {V2[k]}. The weight of each {V2[k]} is
    {W[k] = wx[ix+hx]*wy[iy+hy]} where {k = (ix+hx) + nx*(iy+hy)}.
    
    To obtain one value {V2[k]}, the procedure extracts {ni} sample
    values {y[0..ni-1][k]} and mask values {m[0..ni-1][k]}, where
    {y[i,k]} is the value of image {i} at the point {p[i,k]}, and
    {m[i,k]} is the corresponding inside-outside mask value, for {i}
    in {0..ni-1}. The mismatch {V2[k]} is the variance of the {ni}
    values {y[0..ni-1][k]}, weighted by the corresponding mask values
    {m[0..ni-1][k]}. The variance is estimated assuming that each
    {y[i,k]} has a noise variance {(1 - m[i,k])^2/12}.
    
    This mismatch measure is appropriate when all images have been
    captured using the same instrument with the same light measurement
    scale.It is not adequate, for example, when each image is a
    photograph that has been independently adjusted for brightness and
    constrast, or under different lightings. In these situations,
    one should first apply local or global image normalization filter 
    that eliminates as much of those variations as possible.

    A mismatch function is sensitive only to image /differences/, in the
    sense that it does not change when the same constant {d} gets added
    to all samples of all images. On the other hand, it depends
    quadratically on those differences, so that multiplying all image
    samples by a constant {r} multiplies the result by {r^2}.*/

#endif
