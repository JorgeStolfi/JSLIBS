/* Last edited on 2025-04-27 10:37:37 by stolfi */
/* Reduction of a height difference graph by decimation of low-degree vertices. */ 

#ifndef psgr_shrink_H
#define psgr_shrink_H

/* Created by Rafael F. V. Saracchini */

#include <stdint.h>

#include <psgr_types.h>
#include <psgr.h>

psgr_t *psgr_shrink(psgr_t *gri, int32_t vj_from_vi[], int32_t ej_from_ei[], int32_t indent);
  /* Produces from {gri} a smaller graph {grj} by removing a special set
    of vertices and rearranging their edges so as to maintain
    connectivity and (as far as possible) the same implied height field.
    Like {psgr_shrink_decimate} (q. v.) but operates on a new copy of
    {gri}, thus does not modify {gri}.
    
    Fills {vj_from_vi[0..NVi-1]} and {ej_from_ei[0..NEi-1]} 
    with the mapping of the vertex and edge indices in {gri}, respectively,
    to those in {grj}, as per {psgr_copy} (q. v.).
    
    If {indent} is non-negative, prints diagnostics indented by {indent} columns. */

void psgr_shrink_decimate(psgr_t *gr, int32_t indent);
  /* Deletes from {gr} a subset of the vertices and replaces
    their edges so as to maintain connectivity and (as far as possible)
    the same implied height field.
    
    The vertices and edges are not actually removed from {gr} but only
    disconnected from the graph and marked as deleted. An edge is
    disconnected by removing it from its two initial {onext} rings,
    replacing it as the {.aout} of its endpoints, and marking it
    {psgr_mark_DELETED}. A vertex is removed after all its incident
    edges have been removed, by marking it {psgr_mark_DELETED}.
    
    Upon entry, all vertices and edges must be marked as
    {psgr_mark_CLEAR}.
    
    If {indent} is non-negative, prints diagnostics indented by {indent} columns. */

#endif
