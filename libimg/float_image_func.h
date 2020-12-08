#ifndef float_image_func_H
#define float_image_func_H

/* Procedurally definied float images. */
/* Last edited on 2018-03-04 22:43:30 by stolfilocal */ 

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <sign.h>
#include <float_image.h>

/* 
  A /procedural image/ is a function {func} that returns a float value
  {func(x,y)} for any {double} coordinates {x,y}. 

  When combining procedural images with a pixel array {A} of type
  {float_image_t}, the domain of the latter is assumed to be
  an axis-aligned rectangle {[0 _ NX] × [0 _ NY]} of the plane
  {R^2}, where {NX = A->sz[1]} and {NY = A->sz[2]}.
  
  The pixel in row 0, column 0 of {A} is assumed to be adjacent to the
  origin. More generally, each pixel of {A} with indices {(x,y)} is
  conceptually a square with side 1, whose corners are {(x,y)} and
  {(x+1,y+1)}, and whose center is {(x+0.5,y+0.5)}. */
  
typedef float float_image_func_t(double x, double y);
  /* Type of a function that computes one sample of a procedural image. */
  
/* UTILITY FUNCTIONS */

void float_image_func_get_index_range(double z, double r, int32_t n, int32_t *iLo, int32_t *iHi);
  /* Sets {*iLo} and {*iHi} to the maximum and minimum indices of any
    pixel that may be affected by painting a figure whose real
    coordinates span the range from {z - r} to {z + r}. Takes into
    account that anti-aliasing effectivey enlarges each pixel by 1/2
    pixel all around.
    
    If {n} is not zero, the resulting range is clipped to the range
    {0..n-1}. If {*iLo > *iHi} on return, no pixels are affected. */
  

#endif
