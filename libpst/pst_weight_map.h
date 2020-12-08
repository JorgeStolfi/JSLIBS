#ifndef pst_weight_map_H
#define pst_weight_map_H

/* pst_weight_map.h -- procedures for working with pixel weight maps. */
/* Last edited on 2010-05-04 01:36:51 by stolfi */

#include <bool.h>

#include <float_image.h>

/* WEIGHT MAPS
  
  A /weight map/ is a float-valued image {IW} where each element
  {IW[c,x,y]} is a non-negative weight.
  
  Often the weight {IW[c,x,y]} is the reciprocal of the variance of the 
  noise present in the sample {IM[c',x,y]} of some other channel {c'}
  of some other image {IM}.
*/

float_image_t *pst_weight_map_shrink(float_image_t *IW, bool_t harmonic, int avgWidth);
  /* Given a weight map {IW} for some image {IM}, returns another
    weight map {JW} appropriate for a half-sized version {JM} of the
    image.
    
    Each weight in the reduced map is an average of four weights of
    the given map. The average is the arithmetic mean if {harmonic} is
    false, of the harmonic mean if {harmonic} is true.
    
    If the width of {IW} is odd, the last column is implicitly doubled
    before the image is shrunk. Ditto for the last row, if the height
    is odd. !!! Instead, should assume that pixels outside the domain
    have weight zero. !!!
    
    */
  
float_image_t *pst_weight_map_expand_height_weights(float_image_t *IW);
/*Given a slope weight map {IW}, creates a expanded height weight map where its elements are the 
averaging of equivalent neighbors pixels in SW*/
      
float_image_t *pst_weight_map_slope_to_height(float_image_t *W, bool_t harmonic, int NXV, int NYV);
  /* Given a weight map {W} for a slope map, returns another
    weight map {V} appropriate for the corresponding height map,
    which is assumed to have {NXV} columns and {NYV} rows.
    
    The samples of {W} are assumed to be located at pixel centers,
    that is, sample {W[x,y]} is about the unit square centered at
    the point {(x+0.5,y+0.5)}.  Either the two maps have the same size, or the {V} map has one
    row and one column more than {W}.
    
    If the {V} map is larger, then the samples of {V} are
    associated with vertices of the integer grid.  Namely, 
    each pixel {V[x,y]} of {V} is assumed to be located 
    at the point {(x,y)}. Then sample {V[x,y]} will be the average
    of the {W} samples associated to the four pixels surrounding that point.
    
    If the two maps have the same size, then the samples of {V} are
    assumed to be co-located with the samples of {W}. Then each sample
    {V[x,y]} is the average of nine samples of {W}, associated with
    pixel {[x,y]} and its eight closest neighbors.
    
    The average is the arithmetic mean if {harmonic} is false, of the
    harmonic mean if {harmonic} is true.  Any {W} samples whose
    pixels lie outside {W}'s domain are assumed to be zero. */

/* DEBUGGING */
    
typedef void pst_weight_map_debug_proc_t(int level, float_image_t *W); 
  /* Type of a client-given procedure that may be called
    by recursive integrators to report the weight map used in each scale. */   
float_image_t *pst_weight_map_heights_from_slopes(int NX_Z, int NY_Z,float_image_t *GW);
  
#endif
