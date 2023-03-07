#ifndef frgb_Y_slice_H
#define frgb_Y_slice_H

/* frgb_Y_slice.h - mapping real numbers to paths in the RGB cube */
/* Last edited on 2023-03-06 19:44:32 by stolfi */
/* See the rights and conditions notice at the end of this file. */

#define _GNU_SOURCE
#include <stdint.h>

#include <frgb.h>
#include <bool.h>

void frgb_Y_slice_corners(double z, int32_t *nf_P, frgb_t f[]);
  /* Computes the intersection of the plane {Y=z} with the unit RGB cube,
    where {Y} is the luminance {R*YR + G*YG + B*YB}.  The intersection is a 
    polygon with up to 5 corners. The procedure determines the number {nf}
    of corners, stores the corned into {f[0..nf-1]}, and returns {nf} in {*nf_P}.
    
    If {z < 0} or {z > 1}, the intersection is empty, and {nf} will be zero.
    If {z == 0} or {z == 1}, the intersection is a single point, and {nf}
    will be 1. Otherwise the corners {f[0..nf-1]} will be all distinct,
    and stored in CCW order as seen from the White corner {(1,1,1)}. */
 
#endif

