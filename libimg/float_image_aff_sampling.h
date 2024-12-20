#ifndef float_image_aff_sampling_H
#define float_image_aff_sampling_H

/* Tools for sampling an image with affine deformation. */
/* Last edited on 2024-12-05 10:29:27 by stolfi */ 

#include <bool.h>
#include <r2.h>
#include <i2.h>
#include <r3x3.h>
#include <float_image.h>

r2_t float_image_aff_sampling_choose_step(r3x3_t *M);
  /* Returns the steps in {x} and {y} on {\RR^2} that are adequate to
    sample an image after mapping by a projective map {H} with direct
    matrix {M}. Requires {H} to be affine; that is, the first column of
    {M} should be {[w,0,0]} for some non-zero {w}.
    
    Specifically, determines an orthogonal grid of points on the plane {\RR^2},
    including the origin, that is sufficiently dense for its image under {H} 
    to fully sample an image with sub-pixel accuracy. */
    
i2_t float_image_aff_sampling_grid_size(r2_t dp, double R);
  /* Returns a pair {(NX,NY)} of odd integers that define the size of
    the origin-symmetric sampling grid needed to cover all points of  {\RR^2} within distance
    {R} from the orign.  
    
    Specifically, the sampling point on colum {ix} and row {iy} of this
    sampling grid, for {ix} in {0..NX-1} and {iy} in {0..NY-1}, is assumed to 
    have coordinates {(dp.c[0]*(ix/(NX-1) - 1/2), dp.c[1]*(iy/(NY-1) - 1/2)}. */

#endif
