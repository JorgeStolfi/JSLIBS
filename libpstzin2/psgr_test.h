/* Last edited on 2025-04-27 02:54:53 by stolfi */
/* Tools for testing the graph functions. */

#ifndef psgr_test_H
#define psgr_test_H

#include <stdint.h>

#include <r2.h>

#include <psgr_types.h>
#include <psgr.h>
#include <psgr_path.h>

psgr_t* psgr_test_make_graph(uint32_t nf, uint32_t nh);
  /* Creates a test graph.  
  
    The graph consists mainly of a rectangular grid of vertices with
    holes, and seven special vertices. There are seven holes, each with
    {nh} rows and {nh} cols. The holes leave a frame {nf} vertices thick
    around the whole grid and between every two adjacent holes. These
    numbers must be positive and {nh} must be odd.
    
    For each {k} in {0..6}, in the middle of hole {k} there is a special
    vertex of degree {k}, connected to {k} vertices along the boundary
    of the hole. */

psgr_mark_t psgr_vertex_mark_throw(bool_t deleted);
psgr_mark_t psgr_edge_mark_throw(bool_t deleted);
  /* These procedures return a random mark value suitable for a vertex
    or an edge, respectively.  The result may be {psgr_mark_DELETED}
    only if {deleted} is true. */ 

#endif


