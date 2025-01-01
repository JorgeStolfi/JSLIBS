#ifndef cpk_graph_H
#define cpk_graph_H

/* Undirected graphs represented by neighbor lists. */
/* Last edited on 2024-12-31 14:43:39 by stolfi */

#include <stdint.h>

#include <bool.h>

#include <cpk_basic.h>

/* Get origin and destination of an edge represented as {ui2_t}: */
#define ORG(e) ((uint32_t)(e).c[0])
#define DST(e) ((uint32_t)(e).c[1])

typedef struct cpk_graph_t 
  { uint32_t nV;     /* Number of vertices. */
    uint32_t nE;     /* Number of (undirected) edges. */
    uint32_t *deg;   /* Vertex degrees. */
    uint32_t *fnb;   /* Start of neighbors. */
    uint32_t *nbr;   /* The neighbors. */
  } cpk_graph_t;
  /* An undirected graph in neighbor list format. 
  
    The vertices are represented by indices {0..nV-1}. The neighbors
    of each vertex {v} are {nbr[k..k+m-1]}, where {k = fnb[v]} is the
    index of the first neighbor in {nbr}, and {m = deg[v]} is the
    degree of {v}. Moreover, {fnb[nV] = 2*nE}, and {deg[v]} is always
    equal to {fnb[v+1]-fnb[v]}.  */

cpk_graph_t cpk_gather_neighbors(uint32_t nV, ui2_vec_t E);
  /* Given the edges {E[0..NE-1]} of a graph on {NV} vertices, identified by
  indices {0..NV-1}, returns the neighbor-list representation of that graph. 

  Assumes that the graph is undirected and loop-free, and that each
  undirecte edge is listed exctly once. That is, if the pair {(u,v)}
  is listed in {E}, then {u != v}, and the pair {(v,u)} is not
  listed. */

bool_t cpk_graph_adjacent(cpk_graph_t *G, uint32_t u, uint32_t v);
  /* TRUE iff {u} is a neighbor of {v}. */

uint32_t cpk_graph_count_common_neighbors(cpk_graph_t *G, uint32_t u, uint32_t v, bool_t temp[]);
  /* Number of vertices that are adjacent to both {u} and {v}.
    The {temp} parameter is a client-supplied working vector, with
    at least {G->nV} elements; it must be all FALSE upon entry,
    and will be all FALSE upon exit. */

#endif
