#ifndef pst_graph_from_maps_H
#define pst_graph_from_maps_H

/* Last edited on 2025-03-15 10:17:59 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdint.h>

#include <float_image.h>

#include <pst_graph.h>

pst_graph_t* pst_graph_create_from_gradient_and_weight_maps(float_image_t* IG, float_image_t* IW);
  /* Create a graph {g} from a gradient image {IG} and weight image
    {IW}. The images must have the same number {NX} of columns and {NY}
    of rows. Image {IG} must have two channels, the {X} and {Y}
    derivatives. Image {IW} must have only one channel.
  
    Each vertex {v} of {g} corresponds to a vertex {(v.x,v.y)} of the
    image domain grid, with {v.x} in {0..NX} and {v.y} in {0..NY}. 
    
    The vertex index {v.id} is a unique combination {v.x} and {v.y}.
    Note that this correspondence may be lost in reduced copies of
    the graph. 
    
    Each image pixel with indices {(ix,iy)} for {ix} in {0..NX-1} and
    {iy} in {0..NY-1}, is the average of the corresponding quantity
    (slope or weight) in the unit square with corners {(ix,iy)} and
    {(ix+1,iy+1)}. */

void pst_graph_restore_vertex_index(int32_t id, int32_t NX, int32_t NY, int32_t *ix, int32_t *iy);
  /* Returns the pixel indices {(ix,iy)} given the vertex id computed by
    {pst_graph_compute_vertex_index} and the original image dimensions
    {NX}x{NY}. if {id is not in {0..NX*NY-1}, returns {-1} in {*ix} and
    {*iy}. */

int32_t pst_graph_compute_vertex_index(int32_t ix, int32_t iy, int32_t NX, int32_t NY);
  /* Given a pair {(ix,iy)} of pixel coordinates, computes a unique id, given the dimensions
    {NX}x{NY} of the original image.  If either {ix} is not in {0..NX-1} or {iy} is
    not in {0..NY-1}, returns {-1}. */

#endif
