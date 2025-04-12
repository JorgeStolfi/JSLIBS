/* Last edited on 2025-03-14 19:23:55 by stolfi */
/* Tools for testing the graph functions. */

#ifndef pst_gr_test_H
#define pst_gr_test_H

#include <stdint.h>

#include <r2.h>
#include <pst_gr.h>
#include <pst_gr_path.h>

pst_gr_t* pst_gr_test_make_graph(uint32_t nf, uint32_t nh);
  /* Creates a test graph.  
  
    The graph consists mainly of a rectangular grid of vertices with
    holes, and seven special vertices. There are seven holes, each with
    {nh} rows and {nh} cols. The holes leave a frame {nf} vertices thick
    around the whole grid and between every two adjacent holes. These
    numbers must be positive and {nh} must be odd.
    
    For each {k} in {0..6}, in the middle of hole {k} there is a special
    vertex of degree {k}, connected to {k} vertices along the boundary
    of the hole. */

pst_gr_path_t pst_gr_test_throw_path(r2_t *p0, uint32_t n, r2_t *p1);
  /* Generates a random path that is supposed to start at {p0}
    and end at {p1}.  The path has {n} internal vertices, which are 
    equally spaced points along the line from {p0} to {p1}
    plus random lateral displacements. */

#endif


