#ifndef neuromat_eeg_image_basis_H
#define neuromat_eeg_image_basis_H

/* NeuroMat tools to produce an interpolation basis in image form. */
/* Last edited on 2021-08-24 08:38:46 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <r2.h>
#include <r3.h>
#include <float_image.h>

#include <neuromat_eeg.h>
#include <neuromat_eeg_geom.h>
#include <neuromat_eeg_func_basis.h>
#include <neuromat_eeg_image.h>
    
float_image_t **neuromat_eeg_image_basis_make
  ( int32_t ne, 
    neuromat_eeg_func_basis_eval_t eval, 
    float_image_t *msk, 
    int32_t msub,
    r3_t *srad, 
    r2_t *ictr, 
    r2_t *irad
  );
  /* Computes a list {bas[0..ne-1]} of monochrome float-valued images where 
    {bas[i]} is shows the influence of electrode {i} at every point
    of the scalp.
    
    The values of images {bas[0..ne-1]} at some pixel {p} are proportioal 
    to the values {bval[0..ne-1]} returned by {eval(ne, bval, p3D)} where
    {p3D} is the 3D point of the /idealized scalp/ corresponding to the pixel {p}.

    See {neuromat_eeg_geom.h} for the definition of the /schematic/ and 
    /idealized/ coordinate systems.  Each pixel coordinate {p}
    corresponds to point {p2D = (xs,ys)} of the schematic 2D scalp
    where {xs = (ix + 0.5 - ictr.c[0])/irad.c[0])} and {ys = (iy + 0.5
    - ictr.c[1])/irad.c[1])}.  The schematic point {p2D} is the
    mappe to the /schematic 3D scalp/ (unit sphere) as described in 
    {}, and then to the idealized 3D scalp by stretching each axis {i}
    by the idealized scalp radius {srad.s[i]}.
    
    The image {msk} must be a single-channel mask for the head outline
    (1 inside, 0 outside, possibly antialiased along the border). Each
    basis element will have the same size as {msk} and will be zero
    where {msk} is zero. 
    
    The basis images are antialiased by averaging the basis functions
    evaluated at a grid of {msub} times {msub} points inside each pixel.
    However the edges of the image are NOT antialiased.  The 
    mask {msk} should be used as the opacity channel in the final 
    presentation image. */
    
#endif
