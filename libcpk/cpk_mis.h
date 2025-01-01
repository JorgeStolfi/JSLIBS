#ifndef cpk_mis_H
#define cpk_mis_H

/* Tools for efficient maintenance of independent sets of vertices. */
/* Last edited on 2024-12-31 16:26:23 by stolfi */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <pqueue.h>
#include <bool.h>
#include <vec.h>

#include <cpk_graph.h>

/* 
  MAINTAINING INDEPENDENT SETS
  
  The procedures below manipulate two sets {S} and {Q} of vertices of
  an undirected graph {G} with {nV} vertices. The set {S} is
  /independent/ of itself, namely it contains no pair of adjacent
  vertices. The set {Q} is /{S}-admissible/, namely it is disjoint
  from {S} and contains no neighbor of any vertex in {S}. Note that
  {Q} need not be *all* {S}-admissible vertices of {G}, and {S} need
  not be a maximal independent set in {G-Q}. in particular, {S = Q =
  {}} is a valid pair, as well as {S = {}, Q = {0,..nV-1}}. The two
  sets are representd by a data structure {SQ} of type
  {cpk_mis_state_t}. */

typedef struct cpk_mis_state_t
  { cpk_graph_t *G; /* The graph {G}. */
    uint32_t nS;         /* Cardinality of {S}. */
    uint32_t *S;         /* The elements of {S} are {S[0..nS-1]}, in any order. */
    double *W;           /* {W[v]} is the weight of vertex {v}. */
    double WS;           /* Total weight of {S}. */
    uint32_t *degS;      /* {degS[v]} counts neighbors of {v} that are in {S}. */
    uint32_t *locS;      /* {locS[v]==i} iff {S[i] == v} for some {i} in {0..nS-1}; else {UINT32_MAX}. */
    pqueue_t *Q;         /* Elements of {Q}, sorted by decreasing score. */
    uint32_t *degQ;      /* {degQ[u]} counts neighbors of {u} that are in {Q}. */
  } cpk_mis_state_t;
  
#define locNONE (UINT32_MAX)
  /* If {SQ.locS[v]} is {}, then {v} is not in {SQ.S}. */

/* 
  STATE CREATION AND RECLAMATION */

cpk_mis_state_t *cpk_mis_state_new(cpk_graph_t *G, double W[]);
  /* Allocates a new state {SQ} and its internal data structures,
    large enough for a graph with {nV} vertices. The sets {S} and {Q}
    are empty. Time: O(1), but large.
    
    When reclaiming the record {SQ}, beware that the graph {G} and the
    vector {W} are shared, not copied. */

void cpk_mis_state_free(cpk_mis_state_t *SQ);
  /* Reclaims all internal storage hanging from {SQ} (except the graph
    {G} and the weight vector {W}), and {SQ} itself. Time: O(1). */

/* 
  VERTEX SCORING */

double cpk_mis_score(cpk_mis_state_t *SQ, uint32_t w);
  /* Computes the score of a vertex {w} of {Q}, which determines the
    ordering of {Q} in the priority queue. It is defined as the
    "benefit/cost ratio" of adding {w} to {S}: how much weight we
    would gain (namely, {W[w]}) divided by how many admissible
    vertices we would lose (namely, {1 + degQ[w]}). */

/* 
  STATE INITIALIZATION */

void cpk_mis_reset(cpk_mis_state_t *SQ, bool_t verbose);
  /* Resets {S} and {Q} to empty. Does not free any storage. 
    Time: O(nV). */

void cpk_mis_init(cpk_mis_state_t *SQ, bool_t verbose);
  /* Sets {S} to empty, {Q} to all vertices of {G}. Does not 
    allocate or free any storage. Time: O(nV log(nV)). */
    
/* 
  AUTOMATIC S<->Q MOVES 
  
  The following procedures change {S} by one vertex, and automatically
  change {Q} to match. They are no-ops if the vertex {u} is already in
  the relevant set (for {add}) or outisde it (for {delete}).
  
  These procedures take time O(K*(1+log |Q|)) in the worst case
  where {K} is the sum of the degrees of {u} and its neighbors. */

void cpk_mis_add_indep(cpk_mis_state_t *SQ, uint32_t u, bool_t verbose);
  /* If {u} is not in {S}, appends it to {S}, and adds its weight
    {W[u]} to {WS}. Automatically deletes {u} and all its neighbors
    from {Q}, if they are there. Does not change the order of 
    vertices already in {S}. */

void cpk_mis_delete_indep(cpk_mis_state_t *SQ, uint32_t u, bool_t verbose);
  /* Deletes vertex {u} from {S}, subtracting its weight {W[u]} from
    {WS}. Automatically adds to {Q} the vertex {u} and the neighbors
    of {u}, if they become admissible. If {u} is the very last vertex
    of {S}, the procedure will preserve the order of the remaining
    vertices; otherwise their order may change. */

void cpk_mis_copy_indep(cpk_mis_state_t *SQ, uint32_t nX, uint32_t X[], bool_t verbose);
  /* Substitui o conjunto independente {S} pelo conjunto {X[0..nX-1}},
    atualizando {WS} e as demais estruturas de {SQ}. Automaticametne
    coloca em {Q} o conjunto de todos os vértices admissíveis de {G}. 
    Tempo: O(nV log(nV)). */

/*  
  ONLY-S AND ONLY-Q MOVES
  
  The following procedures update only {S} or only {Q}, and can only
  be used in restricted situations. They are no-ops if the vertex is
  already in the relevant set (for {add}) or outside it (for
  {delete}). */

void cpk_mis_simply_add_adm(cpk_mis_state_t *SQ, uint32_t u);
  /* Adds a vertex {u} to the set {Q}. The vertex must not be in {S}
    nor a neighbor of {S}. Time: O(deg[u]*(1+log |Q|)).  */
 
void cpk_mis_simply_delete_adm(cpk_mis_state_t *SQ, uint32_t u);
  /* Deletes vertex {u} from {Q}. Time: O(deg[u]*(1+log |Q|)). */

void cpk_mis_simply_add_indep(cpk_mis_state_t *SQ, uint32_t u);
  /* Appends vertex {u} to set {S}, and adds its weight {W[u]} to
    {WS}. The vertex must not be in {Q} nor adjacent to any vertex of
    {Q} or {S}. O(deg[u]). */

void cpk_mis_simply_delete_indep(cpk_mis_state_t *SQ, uint32_t u);
  /* Deletes vertex {u} from {S}, and subtracts its weight {W[u]} from
    {WS}. If {u} is not the last vertex of {S}, may change the order
    of the remaining vertices. Time: O(deg[u]). */


#endif
