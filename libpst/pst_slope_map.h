#ifndef pst_slope_map_H
#define pst_slope_map_H

/* pst_slope_map.h -- procedures for working with slope maps. */
/* Last edited on 2025-01-07 12:54:30 by stolfi */

#include <bool.h>
#include <r2.h>

#include <float_image.h>
#include <pst_weight_map.h>
#include <pst_height_map.h>

/* SLOPE MAPS
  
  A /slope map/ is a two-channel float-valued image that gives the gradient
  of some height field {Z(X,Y)} --- that is, its derivatives
  {dZ/dX} and {dZ/dY}.
  
  Specifically, samples {G[0,x,y]} and {G[1,x,y]} of a slope map {G} should be the 
  average value of {dZ/dX} and {dZ/dY}, respectively, inside the square on the {XY} plane 
  with corners {(x,y)} and {(x+1,y+1)}.
  
  !!! TO DO !!!
  
  !!! Make {pst_slope_map_build_integration_system_1} accept a weight map, too.
  
  !!! Generalize to use a 2-channel weight map, with separate
  reliability weights for the X and Y component of the gradient. 
  
*/

r2_t pst_slope_map_get_pixel(float_image_t *G, int32_t x, int32_t y);
  /* Extracts the derivatives {dZ/dX} and {dZ/dY} from channels 0 and
    1 of the pixel in column {x}, row {y} of {G}, and returns them
    as a {r2_t} gradient vector. */

void pst_slope_map_set_pixel(float_image_t *G, int32_t x, int32_t y, r2_t *grd);
  /* Stores the derivatives {dZ/dX} and {dZ/dY}, taken from the
    gradient vector {grd} into channels 0 and 1 of the pixel in column
    {x}, row {y} of {G}. */
    
/* ESTIMATING HEIGHT DIFFERENCES */

/* The following procedures estimate differences of the underlying 
  height field {Z(X,Y)} between
  two points of the integer grid, which are corners of the pixels of 
  the given slope map {G}.   They also return the weight (reliability) 
  of that difference, estimated by interpolating a weight map {W}.
  
  The weight map {W} must have a single channel and the same col and row
  counts as {G}. The value of {W[0,x,y]} is assumed to be the reliability 
  weight of the gradient value {G[c,x,y]}.  The map {W} may be
  {NULL}, in which case the weight is assumed to be 1 in every pixel. */

void  pst_slope_map_get_axial_edge_data
  ( float_image_t* G,
    float_image_t* W,
    int32_t x, int32_t y,
    int32_t axis, int32_t dir,
    double *dP, double *wP
  );
  /* Sets {*d} to the height difference {d} between te height values at
    two adjacent grid points, estimated by interpolating the gradient map {G}.
    Also sets {*wP} to the weight {w} of that estimate.
    
    One of the two pixels in question is {x,y}; the other pixel {x',y'}
    is determined by {axis}, which must be 0 or 1, and {dir} must be
    {+1} or {-1}. Specifically, {x',y'} is {x+dir,y} if {axis} is 0, and
    {x,y+dir} if {axis} is 1.  
    
    The values of {d} and {w} are estimated by interpolating two or
    more pixels of {G} and {W} that are adjacent to the grid edge 
    from {x,y} to {x',y'}. */

void pst_slope_map_get_edge_data
  ( float_image_t* G,
    float_image_t* W,
    int32_t x, int32_t y,
    int32_t ux, int32_t uy,
    double *dP, double *wP
  );
  /* Sets {*d} to the height difference {d} between the height map
    pixels {x,y} and {x',y' = x+ux,y+uy}, estimated by interpolating the gradient map {G}. Also sets
    {*wP} to the weight (reliability) of that difference, estimated by
    interpolating the weight map {W}.
    
    The parameters {ux,uy} must be  {+1}, 0, or {-1}, and must not be both zero.
    If one of them is zero, the values of {d} and {w} are estimated using
    {pst_slope_map_get_axial_edge_data}.  Otherwise the procedure 
    uses {pst_slope_map_get_axial_edge_data} on the four edges of 
    the grid cell with corners {x,y} and {x',y'},
    and combines them appropriately to obtain the estimates of {d} and {w}. */

/* MAP SHRINKING */
 
void pst_slope_and_weight_map_shrink
  ( float_image_t *IG, 
    float_image_t *IW, 
    float_image_t **SG, 
    float_image_t **SW
  );
  /* Given a slope map {IG}, containing the derivative of a height
    function {Z} along X and Y axes, returns another slope map {SG},
    with half the size as {IG}, containing the derivatives of a
    version {SZ} of {IZ} with both dimensions and heights scaled by
    half.  If the given image has size {NX} by {NY},
    the result has size {NX/2} by {NY/2}, rounded up.
    
    A pixel in column {x} and row {y} of the result is conceptually
    centered at the vertex point {(2x+1,2y+1)} of {IG}'s domain.
    
    If {IW} is not null, it should be a monochromatic image with the
    same size as {IG}. Each element {IW[0,x,y]} (which must be
    non-negative) is interpreted as the relative reliability of the
    slopes {IG[c,x,y]}, and is used as a weight in the averaging
    of the slopes. In particular, pixels of {IG} that have zero reliability in {IW} are ignored in
    the averaging.  If {IW} is not null, the preocedure returns also a 
    reduced version {SW} of {IW}.
    
    The images {SG} and {SW} are allocated by the procedure. */

#endif
