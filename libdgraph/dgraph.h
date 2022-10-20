#ifndef dgraph_H
#define dgraph_H

/* Directed or undirected graphs represented as sparse matrices of booleans. */
/* Last edited on 2022-10-20 06:14:19 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <bool.h>
#include <spmat.h>
#include <spmat_io.h>
  
/* GRAPHS */

spmat_typedef(dgraph_t,dgraph,bool_t);
spmat_io_def(dgraph_t,dgraph,bool_t);
  /* A bipartite directed graph represented as a sparse matrix of
    booleans. 
    
    The graph connects two vertex sets, the /row/ or /source/ vertices
    {G.S} with indices {0..G->rows-1} and the /column/ or
    /destination/ vertices {G.D} with indices {0..G->cols-1}. The
    element in row {u} and column {v} is true iff there is a directed
    edge from vertex {u} in {G.S} to the vertex {v} in {G.D}.
    
    The {FALSE} entries are considered trivial and usually not stored.
    The order of the edges in the list {G.e} is not fixed by this
    package, but may be altered by some routines.

    To represent a general directed graph on a set {V} of vertices,
    use {G.rows == G.cols == #V} and implicitly identify sources
    and destinations with the same indices.
    
    To represent an undirected graph, either (1) use a symmetric 
    directed graph, or (2) use a directed graph and implicitly 
    combine entries {(u,v)} and {(v,u)} by some symmetric
    function, e.g. {OR}; or (3) use a directed graph where
    {u <= v} for every edge {(u,v)}, and implicitly
    replace queries on {(u,v)} with {(u > v)} by queries
    on {(v,u)}. */

typedef spmat_size_t dgraph_vertex_count_t;
  /* A count of vertices in a graph. */

typedef spmat_index_t dgraph_vertex_index_t;
  /* The vertices of a graph with {nv} vertices
     are indexed with the integers {0..nv-1}. */
     
#define dgraph_NO_VERTEX_INDEX (spmat_NO_INDEX)
  /* A {dgraph_vertex_index_t} value that means `no such vertex'. */
  
/* EDGES */
  
typedef spmat_count_t dgraph_edge_count_t;
  /* A count of edges in a graph. */

typedef spmat_pos_t dgraph_edge_index_t;
  /* The edges of a graph with {ne} edges
    are indexed with integers {0..ne-1}. */
 
#define dgraph_NO_EDGE_INDEX (spmat_NO_POS)
  /* A {dgraph_edge_index_t} value that means `no such edge'. */
 
/* PROCEDURES SPECIFIC TO GRAPHS: */

dgraph_vertex_count_t *dgraph_degrees(dgraph_t *G, int32_t which);
  /* Compute the degrees of all vertices of {G}, returns a vector
    {deg[0..nv-1]}. If {which} is 0, the vector has {nv = G.rows}
    elements, and {deg[u]} is the out-degree of source vertex {u}.
    If {which} is 1, the vector has {nv = G.cols} elements, and
    {deg[v]} is the in-degree of destination vertex {v}. */
  
dgraph_edge_count_t dgraph_add_undirected_edge
  ( dgraph_t *G, 
    dgraph_vertex_index_t u, 
    dgraph_vertex_index_t v, 
    dgraph_edge_count_t ne
  );
  /* Appends to {G} an edge from {u} to {v}, and, if {u != v}, also an
    edge from{v} to {u}. Assumes that the edges of {G} that have been
    defined are {G.e[0..ne-1]}; adds the new edge(s) starting at
    {G.e[ne]}, and returns the count {ne} incremented by the number
    of edges added. */
    
dgraph_vertex_index_t *dgraph_find_spanning_forest(dgraph_t *G);
  /* The set of components of {G} is represented by an array
    {parent[]} that encodes a maximal-size spanning forest of {G}. For
    each vertex {v} of {G}, {parent[v]} is the parent of {v} in that
    spanning forest; or {v} itself if {v} is a root of that forest. */

dgraph_vertex_count_t *dgraph_count_components_by_size(dgraph_t *G);
  /* Returns a vector {ct} such that {ct[k]} is the
    number of connected components of {G} with exactly {k}
    vertices, for each {k} in {0..G.nv}. Note that {ct}
    has {G.nv+1} elements, and {ct[0]} is always 0. */

dgraph_vertex_index_t dgraph_find_root(dgraph_vertex_index_t parent[], dgraph_vertex_index_t u);
  /* Assuming that {parent[]} is a spanning forest of a graph
    represented as above, returns the root {r} of the tree that
    contains {u}. As a side effect, sets {parent[v] = r} for every
    vertex {v} that is an ancestor of {u} in the current forest,
    including {u}. */

/* LIMITS */

#define dgraph_MAX_VERTEX_COUNT (spmat_MAX_SIZE)
  /* The maximum number of vertices in a graph. */

#define dgraph_MAX_VERTEX_INDEX (spmat_MAX_INDEX)
  /* The maximum valid {dgraph_vertex_index_t} value in any graph. */

#define dgraph_MAX_EDGE_COUNT (spmat_MAX_COUNT)
  /* The maximum number of edges in a graph. */

#define dgraph_MAX_EDGE_INDEX (spmat_MAX_POS)
  /* The maximum valid {dgraph_edge_index_t} value in any graph. */

#endif
