/* Last edited on 2025-01-13 08:11:18 by stolfi */
/* Created by Rafael F. V. Saracchini */

#ifndef pst_img_graph_shrink_H
#define pst_img_graph_shrink_H

/* Graph reduction by decimation of low-degree vertices. */ 

#include <stdint.h>

#include <pst_img_graph.h>
#include <pst_img_graph_integrate.h>

pst_img_graph_t *pst_img_graph_shrink(pst_img_graph_t *ig, int32_t jv_from_iv);
  /* Produces from {ig} a smaller graph {jg} 
    by removing a special set of vertices and rearranging
    their edges so as to maintain connectivity and (as far as possible)
    the same implied height field. 
    
    Vertices and edges of {jg} are renumbered consecutively form 0.
    The procedure fills the array {jv_from_iv[0..ig.NV-1]} with 
    a table such that {jv_from_iv[iv]} is the index in {ig}
    of the vertex whose index in {ig} is {iv}; or {-1},
    if {iv} is one of the deleted vertices. */

#endif
