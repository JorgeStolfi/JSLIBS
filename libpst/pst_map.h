#ifndef pst_map_H
#define pst_map_H

/* pst_map.h -- basic ops for photostereo map images. */
/* Last edited on 2025-03-15 21:27:25 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include <vec.h>
#include <float_image.h>

vec_typedef(pst_map_vec_t,pst_map_vec,float_image_t*);
  /* Defines the type {pst_map_vec_t} as a vector of {float_image_t*} elems. */

void pst_map_ensure_pixel_consistency(float_image_t *A, int32_t wch);
  /* Applies {pst_ensure_pixel_consistency(NC,wch,FALSE,v)} to every
    pixel of the map {A}; where {NC} is the number of channels, and
    {v[0..NC-1]} are the samples of that pixel.
    
    Namely, for exery pixel {[X,Y]}, if {A[wch,X,Y]} exists but is zero,
    or any of the samples {A[c,X,Y]} is not finite, seta {A[wch,X,Y]}
    (if it exists) to zero and all other samples {A[c,X,Y]} to {NAN}.
    Fails if any weight {A[wch,X,Y]} exists but is negative or not
    finite. */
   
void pst_map_interpolate_samples
  ( float_image_t *A,
    int32_t c,
    int32_t wch,
    int32_t x0, int32_t y0,
    int32_t x1, int32_t y1,
    bool_t extrapolate,
    double *vR_P, double *wR_P
  );
  /* Estimates the value {vR} of channel {c} of image {A}
    halfway between the centers of the pixels with indices {x0,y0} and
    {x1,y1}, and its reliability weight {wR}. The pixels must be
    adjacent and consecutive, in that order; either vertically ({x1=x0},
    {y1=y0+1}) or horizontally ({x1=x0+1}, {y1=y0}), but need not be
    within the domain of {A}. The results are returned in {*vR_P} and
    {*wR_P}.
    
    Assumes that channel {wch} of {A}, if it exists, has the reliability
    weighs for the slope values in channel {c}. This weight must be
    a finite non-negative number. If a pixel lies outside the domain of
    the map, its weight is taken to be zero.
    
    The value {vR} is obtained by interpolation or extrapolation of a
    list of up to four samples taken from the same row or column that
    contans both pixels {[x0,y0]} and {[x1,y1]}, as symmetrical
    as possible about the edge {e} separating the
    two given pixels {[x0,y0]} and {[x1,y1]}.
    
    Specifically, the procedure uses samples {v[k]=A[c,xs[k],ys[k]]}
    where {k} ranges in {0..n-1} for some {n} in {0..4}, as well as the
    corresponding weights {w[0..n-1]}. These samples must all have have
    nonzero weight and, if {n} is positive, will include at least one of
    the two given pixels. If possible, the list will have two samples on
    each side of {e}; otherwise it will have one sample on one side of
    {e} and 1 to 3 on the other side; otherwise, if extrapolate is true,
    it may have 3 or 4 samples on only one side of {e}. The number {n} will be
    the largest possible satisfying these conditions and priorities.
    
    In particular, {n} will be zero if both given pixels have zero
    weight, or only one has positive weight but {extrapolate} is false.
    
    In any case, if the best {n} is zero, this function returns {vR=NAN}
    and {wR=0}. Otherwise it procedure computes {vR} and {wR} by calling
    {pst_interpolate_values} (q. v.) with arguments
    {n,v[0..n-1],w[0..n-1],m}. */

#endif
