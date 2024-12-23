#ifndef pst_graph_H
#define pst_graph_H

/* Last edited on 2024-12-22 21:37:34 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <float_image.h>

typedef struct pst_edge_t
  { uint32_t u,v;   /* Indices of origin and destination vertices. */
    double d;       /* Height Difference d of the edge */
    double w;       /* Weight w of the edge */
    int32_t axis;   /* Either 0 for horizontal, 1 for vertical, or {-1} for undefined.  */
  } pst_edge_t;
  /* Represents an edge from vertex {u} to vertex {v}. 
    The {axis} is relevant only while the vertices are pixels
    and every edge is either horizontal or vertical, between adjacent pixels. It is
    {-1} in reduced graphs. */
    
typedef struct pst_vertex_t
  { uint32_t id;                    /* Stores unique id number for each vertex. */
    int32_t x, y;                   /* Pixel indices, or {-1}. */
    int32_t vtxm, vtxp, vtym, vtyp; /* Indices of edges out of vertex, or {-1}. */
  } pst_vertex_t;
 /* Represents a vertex of the graph. 
   The pixel indices {x,y} and the cardinal neighbors {vtxm,vtxp,xtym,vtyp} 
   correspond to the image pixels and pixel adjacencies when the graph is 
   built form an image.  In the reduced graph they are nominal. */
  
struct pst_graph_t
  { pst_vertex_t* vertices;  /* Array that contain the vertices of the graph*/
    pst_edge_t* edges;       /* Array that contain the edges of the graph*/
    uint32_t num_vertex;     /* How much vertices exists in g*/ 
    uint32_t num_edge;       /* How much edges exists in g */
    uint32_t max_vertex;     /* Maximum allowed number of vertices*/
    uint32_t max_edge;       /* Maximum number allowed of edges */
  };

typedef struct pst_graph_t pst_graph_t;

pst_graph_t* pst_graph_new(uint32_t max_vertex, uint32_t max_edge);
  /*Alocates a graph structure with enough space for {max_vertex} vertices and {max_edge} edges*/

void pst_graph_add_vertex(pst_graph_t* g, uint32_t id, int32_t x, int32_t y);
  /* Adds the vertex with id {id} and pixel indexes {(x,y)}, which may be {(-1,-1)},
    Initializes the vertex as if it
    was completely disconected from the graph.  Namely, sets {vtxm,vtxp,vtym,vtyp}
    to {-1}, and does not add any edge. */

void pst_graph_add_edge(pst_graph_t* g, uint32_t u, uint32_t v, double d, double w, int32_t axis);
  /* Adds an edge to the graph from {u} to {v} with  height difference {d} and weight {w}. 
    The vertices need not yet have been added yet to the graph.
    Does not update the fields {vtxm,vtxp,vtym,vtyp} of the vertices.  */

int32_t pst_graph_compute_vertex_index(int32_t ix, int32_t iy, int32_t NX, int32_t NY);
  /* Given a pair {(ix,iy)} of pixel coordinates, computes a unique id, given the dimensions
    {NX}x{NY} of the original image.  If either {ix} is not in {0..NX-1} or {iy} is
    not in {0..NY-1}, returns {-1}. */

void pst_graph_restore_vertex_index(int32_t id, uint32_t NX, uint32_t NY, int32_t *ix, int32_t *iy);
  /* Returns the pixel indices {(ix,iy)} given the vertex id computed by {pst_graph_compute_vertex_index}
    and the original image dimensions {NX}x{NY}.  if {id is not in {0..NX*NY-1},
    returns {-1} in {*ix} and {*iy}. */

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

pst_graph_t* pst_graph_create_from_gradient_and_weight_maps(float_image_t* IG, float_image_t* IW);
  /* Create a graph from a gradient image {IG} and weight image {IW}. */

void pst_graph_write(FILE* arq,pst_graph_t* g);
  /* Writes the graph {g} to {arq} in human-readable format. */

pst_graph_t* pst_graph_shrink(pst_graph_t* g);
  /* Createsa graph {h} that is a shrunk version of {g}.  Roughly every
    two vertices of {g} become one vertex of {h}. */

uint32_t* pst_graph_spanning_forest(pst_graph_t* g);
  /* Computes a spanning forest of {g}. Returns an array
    {emark[0..g.num_edge-1]} where {emark[i]} is the index of the tree
    that has edge {i} of the graph, starting from 1 to {num_tree} where
    {num_tree} is the number of components in the graph. Assumes that
    the graph is grid-like so that the neighbors of each vertex {u} are
    {u.vtxm,u.vtxp,u.vtym,u.vtyp} and every edge has {axis} 0 or 1. */

void pst_graph_interpolate_two_samples
  (  float_image_t* I, float_image_t* W,
     uint32_t c,
     int32_t x0, int32_t y0,
     int32_t x1, int32_t y1,
     double *v, double* w
   );

void pst_graph_free(pst_graph_t* g);

#endif
