/* mst.h -- build minimum-cost spanning trees. */
/* Last edited on 2024-11-20 12:38:40 by stolfi */ 

#ifndef mst_H
#define mst_H

#include <stdint.h>

#include <bool.h>

typedef double mst_arc_cost_t(uint32_t u, uint32_t v);
  /* Type of a function that returns the cost of an arc between
    vertices {u} and {v} of a graph. The result must be non-negative
    and not {NAN}. If there is no such arc, it should return
    {+INF}. */

void mst_build_complete(uint32_t n, mst_arc_cost_t *acost, uint32_t P[], double C[], bool_t verbose);
   /* Builds a minimum-cost spanning rooted forest for a
     graph {G} with vertices {0..n-1}.
     
     The vectors {P,C} must have {n} elements. The function {acost} defines the
     cost of the arc between vertices {u,v}.  The arcs of {G} are those pairs of vertices
     whose cost is finite.
     
     The procedure computes a spanning forest {F} in {G} which has the
     maximum number of edges (i.e. one tree in each connected component of {G}) and minimum total cost
     (which is also the spanning forest with minimum lexicographic list of arc costs). The forest has a root 
     vertex in each component of {G}. The forest {F} is described
     by the vectors {P} and {C}. Namely if {P[u] = u} then {u} is a
     root and {C[u]} is infinite, otherwise {P[u]} is the successor of {u} in the path
     from {u} to the root of its tree and {C[u]} is the (finite) cost of the
     arc from {u} to {P[u]}. */

#endif
