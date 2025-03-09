/* Last edited on 2025-02-11 10:10:28 by stolfi */
/* Created by Rafael F. V. Saracchini */

#ifndef pst_img_graph_from_maps_H
#define pst_img_graph_from_maps_H

/* Improved version of {pst_graph_from_maps.h} that uses {haf.h} for the topology info. */ 

#include <r2.h>
#include <bool.h>
#include <float_image.h>

#include <pst_imgsys.h>
#include <pst_img_graph.h>

pst_img_graph_t* pst_img_graph_from_gradient_and_weight_maps(float_image_t* IG, bool_t add_diags);
  /* Creates a height difference graph {g} from the gradient map {IG}.  
    
    The gradient map {IG} must have three channels; the value of each
    pixel is assumed to be {(dzdx,dzdy,wt)} where {(dzdx,dzdy)} is the
    average gradient inside that pixel, and {wt} is the reliability
    weight of that data, a number between 0 and 1, interpreted as
    proportional to the reciprocal of the variance of the noise in the
    gradient components. If the third channel is missing, all weights 
    are assumed to be 1.
    
    The graph will have one vertex for each pixel of the implied height map
    {OZ}. Each pixel of {OZ} corresponds to a corner of the gradient map
    {IG}; therefore, the height map has {NX+1} cols and {NY+1} rows, and
    the graph {g} will have {NX+1)*(NY+1) vertices.
    
    The graph will have one edge between any two vertices that are
    adjacent either horizontally or vertically.  The edge parameters
    {e.d} and {e.w} are obtained by interpolating the values of {IG} and
    {IW}.  Edges with zero weight are omitted.  If {add_diags} is 
    true, an attempt is made to add diagonal edges in lieu of omitted edges.
    (!!! This option is currently disabled !!!) */

int32_t pst_img_graph_from_maps_get_vertex_index_from_height_map_indices(int32_t ix, int32_t iy, int32_t NX_Z, int32_t NY_Z);
  /* Assumes that {NX_Z,NY_Z} are the col and row counts of the height
    map (NOT of the gradient map). Returns a vertex {id} given the
    height map pixel indices {ix} and {iy}, which must be in {0..NX_Z-1}
    and {0..NY_Z-1}, respectively. If they are outside those ranges,
    returns {-1}. */

void pst_img_graph_from_maps_get_height_map_indices_from_vertex_index(int32_t xy, int32_t NX_Z, int32_t NY_Z, int32_t *ix, int32_t *iy);
  /* Assumes that {NX_Z,NY_Z} are the col and row counts of the height
    map (NOT of the gradient map) and {xy} is the combined pixel indices
    {ix + iy*NX_Z}. Returns the separate indices in {*ix} and {*iy},
    which will be in {0..NX_Z-1} and {0..NY_Z-1}, respectively. The
    input {xy} must be in {0..NX_Z*NY_Z-1}; if not, returns {-1} in
    {*ix} and {*iy}. */


#endif
