#ifndef pst_map_H
#define pst_map_H

/* pst_map.h -- basic ops for photostereo map images. */
/* Last edited on 2025-03-01 19:54:54 by stolfi */

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

#endif
