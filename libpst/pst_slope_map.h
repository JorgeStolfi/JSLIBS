#ifndef pst_slope_map_H
#define pst_slope_map_H

/* pst_slope_map.h -- procedures for working with slope maps. */
/* Last edited on 2025-04-03 17:27:16 by stolfi */

#include <bool.h>
#include <r2.h>

#include <float_image.h>
#include <pst_height_map.h>

/* SLOPE MAPS
  
  A /slope map/ is a float-valued image with either 2 or 3 channels that
  gives the gradient of some height field {Z(X,Y)} --- that is, its
  derivatives {dZ/dX} and {dZ/dY} --- and the respective reliability
  weight.
  
  Specifically, samples {G[0,x,y]} and {G[1,x,y]} of a slope map {G}
  should be the average value of {dZ/dX} and {dZ/dY}, respectively,
  inside the square on the {XY} plane with corners {(x,y)} and
  {(x+1,y+1)}. 
  
  If the map has 3 channels, the value of {G[2,x,y]} should indicate the reliability
  of those values, in the scale from 0 (the gradient data is
  meaningless) to 1 (the gradient data is as accurate as it can be).
  Often the weight {G[2,x,y]} is the reciprocal of the variance of the 
  noise present in the samples {G[0,x,y]} and {G[1,x,y]}.
  
  If the map has only 2 channels, the reliability weight is assumed to be 1. */

r2_t pst_slope_map_get_gradient(float_image_t *G, int32_t x, int32_t y);
  /* Extracts the derivatives {dZ/dX} and {dZ/dY} from channels 0 and
    1 of the pixel in column {x}, row {y} of {G}, and returns them
    as a {r2_t} gradient vector. */

float pst_slope_map_get_weight(float_image_t *G, int32_t x, int32_t y);
  /* Extracts the weight from channel 2 of the pixel in column
    {x}, row {y} of {G}.  If the image has only 2 channels, returns {1.0}. */

void pst_slope_map_set_gradient(float_image_t *G, int32_t x, int32_t y, r2_t *grd);
  /* Stores the derivatives {dZ/dX} and {dZ/dY}, taken from the
    gradient vector {grd} into channels 0 and 1 of the pixel in column
    {x}, row {y} of {G}. */

void pst_slope_map_set_weight(float_image_t *G, int32_t x, int32_t y, float w);
  /* Stores the weight {w} into channel 2 of the pixel in column
    {x}, row {y} of {G}. The image {G} must have three channels. */
    
/* ESTIMATING HEIGHT DIFFERENCES */

/* The following procedures estimate differences of the underlying
  height field {Z(X,Y)} between two points of the integer grid, which
  are corners of the pixels of the given slope map {G}. They use the
  derivatives {dZ/dX} and {dZ/dY} which are assumed to be stored in
  channels 0 and 1 of {G}. The procedures also return the weight
  (reliability) of that difference, estimated by interpolating the
  weight channel {G[2,*,*]}. */

void  pst_slope_map_get_axial_edge_data
  ( float_image_t* G,
    int32_t x, int32_t y,
    int32_t axis, int32_t dir,
    bool_t extrapolate,
    double *dP, double *wP
  );
  /* Sets {*d} to the height difference {d} between the height values at
    two adjacent grid points, estimated by interpolating the gradient map {G}.
    Also sets {*wP} to the weight {w} of that estimate.
    
    One of the two grid vertices in question is {x,y}; the other grid point {x',y'}
    is determined by {axis}, which must be 0 or 1, and {dir} must be
    {+1} or {-1}. Specifically, {x',y'} is {x+dir,y} if {axis} is 0, and
    {x,y+dir} if {axis} is 1.  
    
    The values of {d} and {w} are estimated by interpolating the
    gradients of the two grid cells {x0,y0} and {x1,y1} of {G} that are
    adjacent to the grid edge from {x,y} to {x',y'}, plus two pixels
    {xm,ym} and {xp,yp} that are once removed from them. See
    {pst_map_interpolate_samples}. */

void pst_slope_map_get_edge_data
  ( float_image_t* G,
    int32_t x, int32_t y,
    int32_t ux, int32_t uy,
    bool_t extrapolate,
    double *dP, double *wP
  );
  /* Sets {*d} to the height difference {d} between the height map
    pixels {x,y} and {x',y' = x+ux,y+uy}, estimated by interpolating the
    gradient map {G}. Also sets {*wP} to the weight (reliability) of
    that difference.
    
    The parameters {ux,uy} must be  {+1}, 0, or {-1}, and must not be both zero.
    If one of them is zero, the values of {d} and {w} are estimated using
    {pst_slope_map_get_axial_edge_data}.  Otherwise the procedure
    uses the gradient in the cell with corners {(x,y)} and {(x',y')}. */

/* MAP SHRINKING */
 
float_image_t *pst_slope_map_shrink(float_image_t *IG, double scale);
  /* Given a slope map {IG}, containing the derivatives of a height
    function {Z} along X and Y axes, returns another slope map {JG},
    with half the size as {IG} and same number of channels, 
    containing the derivatives of {Z} sampled at half the 
    resoluton.  The gradient (channels 0 and 1) of the
    result are multipled by the given {scale}.
    
    If the map {IG} has size {NXI} by {NYI}, the result {JG} has size
    {NXJ=NXI/2} by {NYJ=NYI/2}, rounded UP. The number of channels of
    {JG} will be the same as that of {IG} (which must be 2 or 3).
    
    The {scale} factor is typically 1.0.  Since the domain 
    size is reduced by half, the {Z} function values implied by map {JG}
    will be scaled by 1/2.
    
    The reduction is performed by
    {pst_cell_map_shrink(IG,2,NXJ,NYJ,scale)} (q.v.). */

/* ADDING NOISE */

float_image_t* pst_slope_map_merge(float_image_t *GA, float_image_t *GB);
  /* Creates a slope map {G} by combining two slope maps {GA,GB} -- for
    example, one computed numerically, one computed analytically. 
    
    The two maps must have the same size, and either 2 or 3 channels. If
    a map has only two channels, it is assumed to have a third channel
    with all 1.0 reliability weights.
    
    The gradient values (channels 0 and 1) of {G} are the average of the
    corresponding values of {GA,GB}, weighted by the respective
    reliability weights. The reliability weight (channel 2) {w} of {G}
    is the harmonic mean of the weights {wA} and {wB} of {GA} and {GB}.
    In particular, if the result is zero (that is, if {wA} and/or {wB}
    are zero), the weight {w} is zero, and the gradient is set to
    {(NAN,NAN)}. */

/* ADDING NOISE */

void pst_slope_map_perturb(float_image_t *G, double sigma, uint32_t seed);
  /* Adds to each sample of {G} a random value
    with independent Gaussian distribution, mean 0 and
    deviation {sigma}.  If {seed} is nonzero,
    resets the random generator with that seed. */

#endif
