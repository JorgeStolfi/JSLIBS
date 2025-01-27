#ifndef pst_slope_map_H
#define pst_slope_map_H

/* pst_slope_map.h -- procedures for working with slope maps. */
/* Last edited on 2025-01-23 13:20:16 by stolfi */

#include <bool.h>
#include <r2.h>

#include <float_image.h>
#include <pst_weight_map.h>
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

void pst_slope_map_interpolate_two_samples
  ( float_image_t *G,
    int32_t c,
    int32_t x0, int32_t y0,
    int32_t x1, int32_t y1,
    double *vRP, double *wRP
  );
  /* Estimates the value value {vR} of channel {c} (0 or 1) of map {G}
    halfway between the centers of the pixels with indices {x0,y0} and
    {x1,y1}, and its reliability weight {wR}. Returns the results in
    {*vRP} and {*wRP}.
    
    Assumes that channel 2 of {G}, if ot exists,
    gives the reliability of channels 0 and 1 of {G}. */
   
void pst_slope_map_interpolate_four_samples
  ( float_image_t *G,
    int32_t c,
    int32_t x0, int32_t y0,
    int32_t x1, int32_t y1,
    double *vRP, double *wRP
  );
  /* Estimates the value {vR} of channel {c} (0 or 1) of image {G}
    halfway between the centers of the pixels with indices {x0,y0}
    and {x1,y1}, and its reliability weight {wR}.  The pixels must be 
    adjacent either vertically or horizontally.  Returns the 
    results in {*vRP} and {*wRP}. 
    
    Uses two other samples of {G} with indices {xm,ym} and {xp,yp} that
    are collinear with those two points and equally spaced, on both
    sides.  Adjusts the interpolation formula appropritely if
    those samples do not exist,
    
    Assumes that channel 2 of {G} has the reliability weighs for the 
    slope values in channels 0 and 1.  If {G} is null, assumes it is 
    all zeros. */
    
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
    double *dP, double *wP
  );
  /* Sets {*d} to the height difference {d} between the height values at
    two adjacent grid points, estimated by interpolating the gradient map {G}.
    Also sets {*wP} to the weight {w} of that estimate.
    
    One of the two pixels in question is {x,y}; the other pixel {x',y'}
    is determined by {axis}, which must be 0 or 1, and {dir} must be
    {+1} or {-1}. Specifically, {x',y'} is {x+dir,y} if {axis} is 0, and
    {x,y+dir} if {axis} is 1.  
    
    The values of {d} and {w} are estimated by interpolating two or
    more pixels of {G} that are adjacent to the grid edge 
    from {x,y} to {x',y'}. */

void pst_slope_map_get_edge_data
  ( float_image_t* G,
    int32_t x, int32_t y,
    int32_t ux, int32_t uy,
    double *dP, double *wP
  );
  /* Sets {*d} to the height difference {d} between the height map
    pixels {x,y} and {x',y' = x+ux,y+uy}, estimated by interpolating the
    gradient map {G}. Also sets {*wP} to the weight (reliability) of
    that difference.
    
    The parameters {ux,uy} must be  {+1}, 0, or {-1}, and must not be both zero.
    If one of them is zero, the values of {d} and {w} are estimated using
    {pst_slope_map_get_axial_edge_data}.  Otherwise the procedure 
    uses {pst_slope_map_get_axial_edge_data} on the four edges of 
    the grid cell with corners {x,y} and {x',y'},
    and combines them appropriately to obtain the estimates of {d} and {w}. */

/* MAP SHRINKING */
 
float_image_t *pst_slope_map_shrink(float_image_t *IG);
  /* Given a slope map {IG}, containing the derivative of a height
    function {Z} along X and Y axes, returns another slope map {JG},
    with half the size as {IG}, containing the derivatives of a
    version {SZ} of {IZ} with both dimensions and heights scaled by
    half.  If the given image has size {NX} by {NY},
    the result has size {NX/2} by {NY/2}, rounded UP.
    The number of channels will be the same as that of {IG} (2 or 3). 
    
    In fractional index terms, point {(x,y)} of {IG}'s domain gets mapped
    to point {(x/2,y/2)} of {JG}'s domain. Thus the pixel {JG[x,y]} in
    column {x} and row {y} of {JG} conceptually corresponds to the {2x2}
    block of pixels {IG[x',y']} of {IG}, where {x'=2*x+dx}, {y'=2*y+dy},
    and {dx,dy} range in {0..1}.
    
    Thus the value of {JG[c,x,y]}, for {c} in {0..1}, will be the
    the weighted average of the four samples {IG[c,x1,y']} above.  The
    averaging will use the weights of those samples in channel 2 of {IG}.
    If the pixel does not exist in {IG}, its weight is taken to be zero.
    
    If {IG} has three channels, the value of {JG[2,x,y]} will be the
    smallest of the weights of those four samples.  This criterion is
    justified by the assumption that a weight {IG[2,x',y']=0} indicates
    that the height {Z} may be discontinuous in that pixel, in which
    case it will be discontinuous in pixel {x,y} of {JG}. */

/* ADDING NOISE */

void pst_slope_map_perturb(float_image_t *G, double sigma, uint32_t seed);
  /* Adds to each sample of {G} a random value
    with independent Gaussian distribution, mean 0 and
    deviation {sigma}.  If {seed} is nonzero,
    resets the random generator with that seed. */

#endif
