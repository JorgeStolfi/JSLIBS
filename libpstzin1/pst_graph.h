#ifndef pst_graph_H
#define pst_graph_H

/* Last edited on 2024-12-23 14:53:33 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdint.h>
#include <stdio.h>

#include <bool.h>
    
typedef struct pst_vertex_t
  { uint32_t id;                    /* Stores unique id number for each vertex. */
    int32_t x, y;                   /* Pixel indices, or {-1}. */
    int32_t vtxm, vtxp, vtym, vtyp; /* Indices of edges out of vertex, or {-1}. */
  } pst_vertex_t;
 /* Represents a vertex of the graph. The vertex is point {(x,y)} of the
   integer grid {\RZ\times\RZ}. It is incident to at most four edges
   that connect it to adjacent pixels. Specifically, the edges with
   indices {vtxm}, {vtxp}, {vtym}, and {vtyp} connect it to vertices
   {(x-1,y)}, {(x+1,y)}, {(x,y-1)} and {(x,y+1)}, respectively. Each of
   these fields is {-1} is the corresponding edge does not exist. */

typedef struct pst_edge_t
  { uint32_t u, v;  /* Indices of origin and destination vertices. */
    double d;       /* Height Difference of the edge. */
    double w;       /* Weight of the edge. */
    int32_t axis;   /* Either 0 for horizontal, 1 for vertical, or {-1} for undefined.  */
  } pst_edge_t;
  /* Represents an edge from vertex {u} to vertex {v}.
    The {axis} specifies the edge direction. The weight {w} must 
    be positive: zero-weight edges are not stored. */
  
struct pst_graph_t
  { pst_vertex_t* vertices;  /* The vertices of the graph. */
    pst_edge_t* edges;       /* The edges of the graph. */
    uint32_t NV;             /* Vertex count. */ 
    uint32_t NE;             /* Edge count. */
    uint32_t NV_max;         /* Total allocated size of the {vertices} array. */
    uint32_t NE_max;         /* Total allocated size of the {edges} array. */
  };
  /* The vertices are {vertices[0..NV-1]} and the edges are {edges[0..NE-1]}.
    The vertex identifiers are all distinct, but otherwise arbitrary, and may not
    be a consecutive. */

typedef struct pst_graph_t pst_graph_t;

pst_graph_t *pst_graph_new(uint32_t NV_max, uint32_t NE_max);
  /*Alocates a graph structure with enough space for {NV_max}
    vertices and {NE_max} edges*/

void pst_graph_free(pst_graph_t* g);

void pst_graph_add_vertex
  ( pst_graph_t* g,
    uint32_t id,
    int32_t x,
    int32_t y
  );
  /* Adds to {g} the vertex with ID {id} and pixel indexes {(x,y)},
    which may be {(-1,-1)}, Initializes the vertex as if it was
    completely disconected from the graph. Namely, sets
    {vtxm,vtxp,vtym,vtyp} to {-1}, and does not add any edge. */

void pst_graph_add_edge(pst_graph_t* g, uint32_t u, uint32_t v, double d, double w, int32_t axis);
  /* Adds an edge to the graph from {u} to {v} with  height difference {d} and weight {w}. 
    The vertices need not yet have been added yet to the graph.
    Does not update the fields {vtxm,vtxp,vtym,vtyp} of the vertices.  */

void pst_graph_get_edge_difference_and_weight(pst_graph_t* g, uint32_t v, int32_t e, double*d, double* w);
   /* Given a vertex index {v} and an edge index {e}, returns the height difference 
     {d} and weight {w} when {e} traversed from {v} to the other endpoint.  Namely, if {v} is the origin of the
     edge, returns {e.d}; if {v} is the destination, returns {-e.d}.
     Fails if {v} does not belong to the edge. */

void pst_graph_vertex_get_neighbour
  ( pst_graph_t *g,
    uint32_t vt0,
    int32_t dx, int32_t dy,
    int32_t *vt1,
    double *d,	
    double *w	
  );
  /* Given a vertex index {vt0} and a direction {(dx,dy)} where {dx,dy}
    are in {-1..+1}, returns in {*vt1} the index of the vertex that lies
    in that direction from {vt0}, provided that there is a 1- or 2-edge
    path between the two vertices. In that case, also returns height
    difference {d} and weight {w} of the edge. If there isn't an edge
    between the two vertices, returns {-1} in {*vt1} and {d} and {w}
    values will be 0. Fails if the given direction is bigger than
    allowed or if it is (0,0) */

void pst_graph_update_neighbours(pst_graph_t* g);
  /* Updates the neighbours {vtxm,vtxp,vtym,vtyp} of all vertices,
    according to existing edges.  Only works for the 
    graph is derived from an image, so that each edge connects
    pixels that are adjacent horizontally or vertically. 
    Will fail if this is not the case. */

bool_t pst_graph_check_consistency(FILE* wr,pst_graph_t* g);
  /* Checks the graph consistency using the edge's list. If {wr} is not null,
    a humam-readable report is written to it. */

void pst_graph_write(FILE* wr, pst_graph_t* g);
  /* Writes the graph {g} to {wr} in human-readable format. */

uint32_t* pst_graph_spanning_forest(pst_graph_t* g);
  /* Computes a spanning forest of {g}. Returns an array
    {emark[0..g.NE-1]} where {emark[i]} is the index of the tree
    that has edge {i} of the graph, starting from 1 to {num_tree} where
    {num_tree} is the number of components in the graph. Assumes that
    the graph is grid-like so that the neighbors of each vertex {u} are
    {u.vtxm,u.vtxp,u.vtym,u.vtyp} and every edge has {axis} 0 or 1. */

#endif
