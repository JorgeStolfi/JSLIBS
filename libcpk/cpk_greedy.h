#ifndef cpk_greedy_H
#define cpk_greedy_H

/* Heuristics for finding large independent sets. */
/* Last edited on 2024-12-31 13:23:47 by stolfi */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <bool.h>
#include <pqueue.h>
#include <vec.h>

#include <cpk_graph.h>
#include <cpk_mis.h>

/* 
  FINDING LARGE INDEPENDENT SETS
  
  These procedures take a graph {G} whose vertices {V[0..nV-1]} are
  points of {\Real^2}, with real weights {W[0..nV-1]}; and return an
  independent subset {S} of the vertices (as a lis of vertex indices).
  
  The procedure try to optimize the total weight of {S}, defined 
  as {WS = SUM { W[S[i]] : i = 0..S.ne}}.
  
  The parameter {seed} is used if needed to initialize random number
  gemerators or to select arbitrary internal parameters. */

uint32_vec_t cpk_greedy_find_indep_set
  ( cpk_graph_t *G,  /* The incompatibilty graph. */
    double W[],      /* Weight of each vertex. */
    bool_t verbose   /* TRUE prints debugging diagnostics. */
  );
  /* Uses a greedy heuristic: at each step, the next
    vertex {v} added to {S} is the one that is still a candidate
    (meaning that it is compatible with all previous vertices) and
    maximizes the ratio {W[v]/(1+K[v])} where {K[v]} is the number of
    neighbors of {v} that are still candidates. */
/*  
  PARTIAL MAXIMIZATION */

void cpk_greedy_maximize(cpk_mis_state_t *SQ, bool_t verbose);
  /* Greedily expands the set {S} (in {SQ}) to a maximal independent
    set of vertices. Namely, repeats takes the maximum scored element
    {u} of {Q} and adds it to {S}, updating {Q} and {WS} as
    appropriate, until {Q} becomes empty. The new vertices, if any,
    are appended to {S}, without disturbing the old ones. */

void cpk_greedy_maximize_random(cpk_mis_state_t *SQ, bool_t verbose);
  /* Similar to {cpk_greedy_maximize}, except that, instead of always 
    choosing the maximum-scored element, selects randomly among 
    all the queue elements, biased towards highest scores. */

void cpk_greedy_maximize_lookahead(cpk_mis_state_t *SQ, bool_t verbose, bool_t temp[]);
  /* Expands {S} to a maximal independent set of vertices, by a greedy
    algorithm with one-vertex lookahead. Namely, at each iteration
    looks for two mutually independent vertices in {Q} whose combined
    weight/cost ratio (accounting for neighborhood overlap) is
    maximum; and then inserts the best vertex {u} of the pair in {S},
    updating {Q} and {WS} as appropriate. This process is repeated until {Q}
    becomes empty. The new vertices, if any, are appended to {S},
    without disturbing the old ones. 
    
    The {temp} parameters is a scratch area with at least {SQ.G.nV} elements
    that must be provided by the client. It must be all FALSE upon entry
    and will be all FALSE on return. */

#endif
