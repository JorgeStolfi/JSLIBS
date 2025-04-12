/* Last edited on 2025-03-11 07:08:11 by stolfi */
/* Reduction of a height difference graph by decimation of low-degree vertices. */ 

#ifndef pst_gr_shrink_H
#define pst_gr_shrink_H

/* Created by Rafael F. V. Saracchini */

#include <stdint.h>

#include <pst_gr.h>
#include <pst_gr_integrate.h>

pst_gr_t *pst_gr_shrink(pst_gr_t *gri, int32_t vj_from_vi);
  /* Produces from {gri} a smaller graph {grj} 
    by removing a special set of vertices and rearranging
    their edges so as to maintain connectivity and (as far as possible)
    the same implied height field. 
    
    Vertices and edges of {grj} are renumbered consecutively form 0.
    The procedure fills the array {vj_from_vi[0..gri.NV-1]} with 
    a table such that {vj_from_vi[vi]} is the index in {gri}
    of the vertex whose index in {gri} is {vi}; or {-1},
    if {vi} is one of the deleted vertices. */

#endif
