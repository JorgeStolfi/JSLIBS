#ifndef float_image_aff_extract_H
#define float_image_aff_extract_H

/* Tools for locally comparing two images with affine deformation. */
/* Last edited on 2020-11-05 23:03:35 by jstolfi */ 

#include <bool.h>
#include <r2.h>
#include <i2.h>
#include <r2_aff_map.h>
#include <float_image.h>

float_image_t *float_image_aff_extract(float_image_t *img, r2_aff_map_t *A, r2_t dp, i2_t size);
  /* Returns the image of a feature as considered by {float_image_aff_compare},
    weighted by the fuzzy circular window. 
    
    Specifically, each pixel of the result is {w(p)*img(A(p))}, where
    {p} is a point of the plane {\RR^2}, and {w(p)} is a two-dim
    Gaussian-bell weight function with mean 0 and deviation 1.
    
    The point {p} is varied over {\RR^2} in an orthogonal grid. The grid
    is symmetric about the origin and has {size.c[j]} samples along each
    axis {j}, with step {dp.c[j]}. The grid sizes must be odd and
    positive. The image pixels are interpolated with a C1 scheme. Thus,
    the pixel on colunm {ix[0]} and row {ix[1]} of the resulting image
    corresponds to the sampling point with coordinate {j} equal to
    {dp.c[j]*(ix[j]/(size.c[j]-1) - 1/2)}, where {ix[j]} is in
    {0..size.c[j]-1}. */

#endif
