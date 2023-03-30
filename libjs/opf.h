/* opf.h -- build forests of minimum-cost paths. */
/* Last edited on 2023-03-18 11:16:38 by stolfi */ 

#ifndef opf_H
#define opf_H

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>

typedef double opf_arc_cost_t(int32_t u, int32_t v);
  /* Type of a function that returns the cost of an arc between vertices {u} and {v} of a graph.
    The result must be non-negative and not {NAN}. If there is no such arc, it should 
    return {+INF}. */

typedef double opf_path_cost_t(double Ca, double Cp);
  /* Type of a function that returns the cost of a path that consists of
     one arc with cost{Ca} followed by a path of cost {Cp} It must never
     be smaller than {Ca} and {Cp}. */

void opf_build_complete(int32_t n, double C[], opf_arc_cost_t *acost, opf_path_cost_t *pcost, int32_t P[], int32_t R[], bool_t verbose);
   /* Builds a spanning rooted optimum-path forest for the complete
     graph {G} with vertices {0..n-1}.
     
     The vectors {C,P,R} must have {n} elements. On input, the values
     of {C[u]} should be cost of the trivial path {(u)}; it must be
     non-negative, possibly {+INF}. The function {acost} defines the
     cost of each arc, and together with {pcost} defines the cost of
     non-trivial paths.
     
     The procedure computes a forest {F} in {G} that, for each vertex
     {u}, contains a path that starts at {u} and has minimum cost
     among all such paths. That path is /the optimum path of {u}/,
     {OP(u)} Each tree in this this forest contains a single root
     vertex {r} such that {OP(r)} is trivial, and which is the
     terminus of {OP(u)} all {u} in the same tree.
     
     The resulting forest is described by the vectors {P}, {C}, and
     {R}. Namely if {P[u] = u} then {u} is a root, otherwise {P[u]} is
     the successor of {u} in the path {OP(u)}. For every {u}, {C[u]}
     will be the (possibly infinite) cost of {OP(u)}, and {R[u]} will
     be its terminus.
     
     Currently uses a brute-force algorithm, very far from the theoretical
     optimum. */

#endif
