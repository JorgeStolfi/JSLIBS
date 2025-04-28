/* Last edited on 2025-04-27 02:57:58 by stolfi */
/* Creation of a height difference graph from a slope (gradient) map. */

#ifndef psgr_from_slope_map_H
#define psgr_from_slope_map_H

/* Created by Rafael F. V. Saracchini */

#include <stdint.h>

#include <bool.h>
#include <float_image.h>

#include <pst_imgsys.h>

#include <psgr_types.h>
#include <psgr.h>

psgr_t* psgr_from_slope_map(float_image_t* IG, bool_t extrapolate, bool_t add_diags);
  /* Creates a height difference graph {gr} from the gradient map {IG}.  
    
    The gradient map {IG} must have three channels; the value of each
    pixel is assumed to be {(dzdx,dzdy,wt)} where {(dzdx,dzdy)} is the
    average gradient inside that pixel, and {wt} is the reliability
    weight of that data, a non-negative number, interpreted as
    proportional to the reciprocal of the variance of the noise in the
    gradient components. If the third channel is missing, all weights 
    are assumed to be 1.
    
    The graph will have one vertex for each pixel {OZ[x,y]} of the
    implied height map {OZ}. Each pixel of {OZ} corresponds to a corner
    of the gradient map {IG}.  Therefore, the height map has {NX+1} cols
    and {NY+1} rows, and the graph {gr} will have {NV=(NX+1)*(NY+1)}
    vertices. The pixel indices {x,y} will be saved as the
    {psgr_vertex_pixel} attribute of the vertex.
    
    The graph will have one edge between any two vertices that are
    adjacent either horizontally or vertically.  The edge parameters
    {e.d} and {e.w} are obtained by interpolating the values of {IG} and
    {IW}.  Edges with zero weight are omitted. 
    
    If {add_diags} is true, an attempt is made to add diagonal edges in
    lieu of omitted edges. (!!! This option is currently disabled
    !!!) */

#endif
