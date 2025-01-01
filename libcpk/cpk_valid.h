#ifndef cpk_validate_H
#define cpk_validate_H

/* Procedures to check the validity of circle packings. */
/* Last edited on 2024-12-31 16:02:29 by stolfi */ 

#include <stdint.h>
#include <stdio.h>

#include <bool.h>
#include <r2.h>

#include <cpk_main.h>
#include <cpk_graph.h>
#include <cpk_basic.h>

/* 
  VALIDATION PROCEDURES
  
  These procedures run consistency tests on the arguments as described below.
  If the tests fail, the procedures print an error message and return FALSE.
  Typical use: {assert(cpk_check_...(...))}. */

bool_t cpk_check_edges (r2_vec_t *V, ui2_vec_t *E, double dMin);
  /* Checks whether the edges in {E} are indeed all the 
    pairs of vertices from {V} that are closer than {dMin}
    to each other.  Time is {O(N^2 + M)} where {N = V.ne},
    {M = E.ne}. */

bool_t cpk_check_indep_set(cpk_graph_t *G, uint32_t nJ, uint32_t *J, bool_t *sel);
  /* Checks whether the integers {J[0..nJ-1]} are all distinct, valid, 
    and non-adjacent vertices of {G}. If {sel} is not NULL,also
    checks whether {sel[v]} is true if and only if {v} ocurrs
    in {J}. */

bool_t cpk_check_solution
  ( cpk_domain_t *C,
    cpk_policy_t *P,
    r2_vec_t *V,
    uint32_vec_t *J
  );
  /* Checks whether {V,J} is a valid solution to the RadCom problem 
    instance described by {C} and {P}.
    
    The vector {V} should contain the EUTM X/Y coordinates of the
    candidate auction points (grid points). The vector {J} should
    contain the indices of the vertices which are being proposed as
    auction centers.
    
    In particular, for each proposed auction {J[i]}, finds the nearest
    neighbor {J[k]} with {j\neq i}, and checks whether it is valid.
    (Presently this is the only check done by the procedure.) */

ui2_vec_t cpk_slow_get_proximity_edges (r2_vec_t *V, double dLim);
  /* Computes the proximity (= incompatibility) edges for {V}, that
    is, all pairs {(i,j)} such that {i<j} and {dist(V[i],V[j]) <
    dLim}. Does it in the slow but sure way; useful for validating
    faster algorithms. */
   
#endif
