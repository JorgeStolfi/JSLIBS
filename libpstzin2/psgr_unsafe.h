/* Last edited on 2025-04-27 02:22:34 by stolfi */
/* Private definitions and unsafe functions for {ps_gr_t} structures. */

#ifndef psgr_unsafe_H
#define psgr_unsafe_H

#include <stdint.h>

#include <r2.h>
#include <i2.h>
#include <psgr_types.h>
#include <psgr_path.h>

/* This interface should not be imported by programs that use this
  library. Calling the functions in this module or directly manipulating
  the fields of {psgr_t}, {psgr_vertex_data_t}, or {psgr_edge_data_t}
  records may easily leave the structure in an inconsistent state. This
  module is provided mostly for internal use of other routines in this
  library package. */

typedef struct psgr_vertex_data_t
  { psgr_arc_t aout;      /* One of the arcs out of the vertex, or {psgr_arc_NONE}. */
    i2_t pixel;           /* Col and row indices into the underlying image. */
    psgr_mark_t vmark;    /* Vertex state mark. */
  } psgr_vertex_data_t;
  /* Attributes of a vertex of a graph. 
  
    The {.pixel} field has the indices of the pixel of an image (such as a
    height map) that are associated to this vertex.
    
    The field {.aout} is one of the arcs that leave the vertex, or is
    {psgr_arc_NONE} if the vertex is isolated.
  
    The {.vmark} is used internally when shrinking or enumerating
    connected parts of the the graph. The value {psgr_mark_CLEAR} is
    normal in most vertices and situations. The other values are used
    temporarily by the graph shrinking algorithm. Vertices marked
    {psgr_mark_DELETED} should be isolated vertices (with degree
    zero) that are about to be eliminated from the graph. This field may
    define the vertex color when plotting. */

typedef struct psgr_edge_data_t
  { psgr_vertex_t org[2];     /* Index of the origin vertex in each direction. */
    psgr_arc_t enext[2];      /* Index of the next arc (oriented edge) around each endpoint. */
    psgr_arc_t eprev[2];      /* Index of the previous arc (oriented edge) around each endpoint. */
    double delta;             /* Height difference in the base direction. */
    double weight;            /* Edge weight (for both directions). */
    psgr_path_t path;         /* Path for plotting the edge in the base direction.  */
    psgr_mark_t emark;        /* Used for connected component enumeration. */
  } psgr_edge_data_t;
  /* Attributes of an (undirected) edge of the graph.
  
    The {.path} field defines a curve for plotting the edge. The curve starts
    at {org[0]}, goes near the vertices of {path}, and ends at {org[1]}.
    
    Let {a} be the base arc of the edge. The field {.enext[0]} is the
    arc that follows {a} in the list of arcs that leave the origin
    vertex {org[0]}, in counterclockwise order. The field {.enext[1]} is
    the arc that follows {psgr_arc_sym(a)} in the list of arcs that leave
    the destination vertex {org[1]}, in counterclockwise order. The
    fields {.eprev[0]} and {.eprev[1]} are the same as {.enext[0]} and
    {.enext[1]}, except that the order is clockwise instead of
    counterclockwise.
    
    The {.emark} field is used internally when scanning the graph. The
    value {psgr_mark_CLEAR} is normal in most edges and situations.
    The other values are used temporarily by the graph shrinking
    algorithm. Edges marked {psgr_mark_DELETED} should be
    disconnected from the graph (not part of any {onext} ring) that are
    about to be eliminated from the graph. This field may defines the
    edge color and style when plotting. */

typedef struct psgr_t
  { uint32_t NV;                  /* Actual number of of vertices. */
    uint32_t NE;                  /* Actual number of (undirected) edges. */
    uint32_t NV_max;              /* Allocated space for vertices */
    uint32_t NE_max;              /* Allocated space for edges */
    psgr_vertex_data_t *vdata;    /* Vertex data, indexed {0..NV-1}*/
    psgr_edge_data_t *edata;      /* Edge data, indexed {0..NE-1} */
    /* Connection to height map: */
    uint32_t NX;                  /* Col count of height map. */
    uint32_t NY;                  /* Row count of height map. */
    psgr_vertex_t *vix;           /* Maps height map indices to vertex indices. */
  } psgr_t;
  /* Private spec of a graph structure header.
  
    The vertex data records of a graph {gr} are {gr.vdata[0..gr.NV-1]},
    and the edge data records are {gr.edata[0..gr.NE-1]}.
    
    The vertices are associated with a subset of the points {(x,y)} of
    an integer grid with {x} in {0..gr.NX-1} and {y} in {0..gr.NY-1}.
    The table table {gr.vix} will have {gr.NX*gr.NY} elements, and
    {gr.vix[gr.NX*y + x]} will be the index of the vertex corresponding
    to pixel {[x,y]}; or {psgr_vertex_NONE} if there is no such
    vertex. */
    
/* WARNING: The procedures in this interface that modify the structure 
  may leave it in an inconsistent state.  */

void psgr_unsafe_arc_set_org(psgr_t *gr, psgr_arc_t a, psgr_vertex_t org);
  /* Sets the vertex with index {org} as the origin of the arc {a}
    (which must not be {NULL}. */
 
void psgr_unsafe_arc_set_onext(psgr_t *gr, psgr_arc_t a, psgr_arc_t b);
  /* Sets the {enext} field of the edge of {a} so that {psgr_arc_onext(a)} becomes {b}. 
    Other {enext} and {oprev} fields must be set too to keep the structure consistent. */
    
void psgr_unsafe_arc_set_oprev(psgr_t *gr, psgr_arc_t a, psgr_arc_t b);
  /* Sets the {enext} field of the edge of {a} so that {psgr_arc_oprev(a)} becomes {b}. 
    Other {enext} and {oprev} fields must be set too to keep the structure consistent. */
    
void psgr_unsafe_vertex_set_out_arc(psgr_t *gr, psgr_vertex_t v, psgr_arc_t a);
  /* Sets {a} as the representative member of the onext ring of {v},
    that will be returned by {psgr_vertex_out_arc}.
    If {a} is not {psgr_NONE}, the origin of {a} must be {v}. */

void psgr_arc_unsafe_split
  ( psgr_t *gr,
    psgr_arc_t a,
    psgr_edge_t *e_P,
    psgr_dir_bit_t *db_P,
    psgr_edge_data_t **ed_P
  );
  /* Given an arc {a}, obtains its underlying undirected edge {e}, its
    direction bit {db}, and (if {gr} is not {NULL}) the address {ed} of
    the corresponding edge data record {gr.edata[e.ix]}. These values
    are returned in {*e_P}, {*db_P}, and {*ed_P}, if they are not
    {NULL}. Also checks that {a} is not {psgr_arc_NONE} and (if {gr}
    is not {NULL}) that {e} is a valid edge index for {gr}.
    
    Beware that the procedures that add new edges to the graph may force a 
    re-allocation of the {gr.edata} table.  In that case, the previously 
    obtained edge data record address {ed} may become invalid, and assignig to it
    may beak the heap. */

psgr_arc_t psgr_unsafe_add_unlinked_edge(psgr_t *gr, psgr_vertex_t org, psgr_vertex_t dst, psgr_path_t P);
  /* Adds to {gr} a new edge {e} connecting existing vertices {org} and
    {dst}. Returns the arc {a} of that edge that is oriented from {org}
    to {dst}.
  
    The new edge index {e.ix} will be {NE-1} where {NE} is
    {psgr_num_edges(gr)} on exit.
    
    The delta and weight of {a} will be set to {NAN}, and the path of
    the arc {a} will be {P}. The edge's mark is set to
    {psgr_mark_CLEAR}.
    
    This procedure does NOT insert {a} and its symmetrical into the
    onext rings of {org} and {dst}, and in fact will not change any
    onext rings. The pointers {psgr_arc_onext(a)},
    {psgr_arc_oprev(a)}, {psgr_arc_dnext(a)}, and
    {psgr_arc_dprev(a)} will all be undefined ({psgr_arc_NONE}). 
    It also does NOT change the {.aout} fields of {org} and {dst},
    even if they are {psgr_arc_NONE}. */

#endif


