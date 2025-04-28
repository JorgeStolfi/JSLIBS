/* Last edited on 2025-04-27 10:44:20 by stolfi */
/* Creating a deep copy of a {psgr_t} graph structure. */ 

#ifndef psgr_copy_H
#define psgr_copy_H

/* Created by Rafael F. V. Saracchini */

#include <stdint.h>

#include <r2.h>

#include <psgr_types.h>
#include <psgr_path.h>
#include <psgr.h>

#include <psgr_copy.h>
  
psgr_t* psgr_copy(psgr_t *gri, int32_t vj_from_vi[], int32_t ej_from_ei[]);
  /* Creates a copy {grj} from {gri}, excluding all vertices and edges
    that are marked as {psgr_mark_DELETED}. 
    
    New copies are made of all records of unmarked vertices and edges,
    including path records and marks. In the new copy, vertices and
    edges are renumbered consecutively starting from 0, and their marks
    are all set to {psgr_mark_CLEAR}.
    
    The procedure fills the array {vj_from_vi[0..NVi-1]}, where {NVi} is
    the vertex count of {gri}, with a table such that
    {vj_from_vi[vi_ix]} is the new index of the vertex whose inde was
    originally {vi_ix}. Similarly, it fills the table
    {ej_from_ei[0..NEi-1]}, where {NEi} is the edge count of {gri} that
    gives the new index of each undeleted edge of {gri} in {grj}. In
    both cases, an entry of {-1} means the vertex or edge of {gri} was
    deleted.
    
    The tables {vj_from_vi} and {ej_from_ei} may be {NULL} if no
    vertices of {gri} are marked for deletion, since in that case the
    vertex and edge indices will be preserved. The procedure will fail
    if either table is {NULL} but {gri} contains some elements of the
    corresponding type that are marked {DELETED}. */

#endif


