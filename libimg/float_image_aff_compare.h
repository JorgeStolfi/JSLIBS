#ifndef float_image_aff_compare_H
#define float_image_aff_compare_H

/* Tools for locally comparing two images with affine deformation. */
/* Last edited on 2022-10-19 08:27:23 by stolfi */ 

#include <bool.h>
#include <r2.h>
#include <i2.h>
#include <hr2.h>
#include <float_image.h>

double float_image_aff_compare
  ( float_image_t *img1,
    hr2_pmap_t *A1,
    float_image_t *img2,
    hr2_pmap_t *A2,
    r2_t *stepP,
    i2_t *sizeP
  );
  /* Computes the local discrepancy squared between two images {img1}
    and {img2} modified by two invertible affine maps {A1} and {A2}.
    
    Specifically, the procedure returns the weighted-average square
    pixel dscrepancy
    
      {INT{p: |f1(p)-f2(p)|^2 w(p)} / INT{p: w(p)}}
    
    where {p} varies over {\RR^2}, {f1(p} = img1(A1(p))}, {f2(p} =
    img2(A2(p))}, and {w(p)} is a two-dim Gaussian-bell weight function
    with mean 0 and deviation 1.
    
    The space {\RR^2} is sampled on an orthogonal grid, with sufficient
    density along each axis to sample both {img1} and {img2} with
    sub-pixel accuracy. The pair of sampling steps is returned in
    {*stepP}, if not NULL. The image pixels are interpolated with a C1
    scheme, so the result is C1 on the coefficients of {A1} and {A2}.
    
    The sampling grid has an odd number of points along each axis, and
    is symmetric about the origin. The grid size is large enough to
    include all samples with significant weight, and is returned in
    {*sizeP}, if not NULL. */

#endif
