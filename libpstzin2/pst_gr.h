/* Last edited on 2025-03-14 20:15:30 by stolfi */
/* Graph encoding of gradient maps for slope-to-height integration. */

#ifndef pst_gr_H
#define pst_gr_H

/* Improved version of {pst_graph.h} that emulates {haf.h} for the topology info. */ 
/* Created by Rafael F. V. Saracchini */

#include <stdint.h>

#include <r2.h>
#include <pst_gr_path.h>

/* This interface defines a data structure to represent a planar graph properly embeded 
  in the plane, where each directed edge specifies the difference between
  the heights of the origina and destination vertices, and a reliability weight 
  for the same. */

typedef uint32_t pst_gr_vertex_t;
  /* Index of a vertex in a graph. */
  
typedef uint32_t pst_gr_edge_t;
  /* Index of an (unoriented) edge in a graph. */
  
typedef uint32_t pst_gr_arc_t;
  /* Index of an arc (oriented edge) in a graph. It is {2*ei + db} where
    {ei} is the index of the unoriented edge, and {db} is a /direction
    bit/ bit that specifies the direction of the arc. The /base arc/ of
    the edge is the one with {db=0}. */
   
typedef uint8_t pst_gr_dir_bit_t;
  /* The direction bit of an arc (oriented edge).  Actually, it can be only 0 or 1. */
   
#define pst_gr_NONE ((pst_gr_arc_t)UINT32_MAX)
  /* A value for a {pst_gr_vertex_t}, {pst_gr_arc_t}, or {pst_gr_edge_t}
    that indicates "no such item". */
   
typedef uint32_t pst_gr_mark_t;
  /* The mark field of an edge or vertex record. */

#define pst_gr_UNMARKED ((pst_gr_mark_t)UINT32_MAX)
  /* A null value for a {pst_gr_mark_t} field. */

typedef struct pst_gr_vertex_data_t
  { pst_gr_arc_t aout;      /* One of the arcs out of the vertex, or {pst_gr_NONE}. */
    int32_t x;              /* Orignal col index in height map, or {-1} if none. */
    int32_t y;              /* Orignal row index in height map, or {-1} if none. */
    r2_t coords;            /* Plot coordinates. */
    pst_gr_mark_t vmark;    /* Used during enumeration of connected components. */
  } pst_gr_vertex_data_t;
  /* Attributes of a vertex of a graph.  The {coords} are used only for plotting the graph.
  
    The fields {x,y}, when not {-1}, are indices of the pixel of a height map that
    are associated to this vertex.
    
    The field {aout} is one of the arcs that leave the vertex, or is {pst_gr_NONE}
    if the vertex is isolated.
  
    The {vmark} is used internally when shrinking or enumerating connected parts
    of the the graph.  It also defines the vertex color when plotting. */

typedef struct pst_gr_edge_data_t
  { pst_gr_vertex_t org[2];     /* Index of the origin vertex in each direction. */
    pst_gr_arc_t enext[2];      /* Index of the next arc (oriented edge) around each endpoint. */
    pst_gr_arc_t eprev[2];      /* Index of the previous arc (oriented edge) around each endpoint. */
    double delta;               /* Height difference in the base direction. */
    double weight;              /* Edge weight (for both directions). */
    pst_gr_path_t path;         /* Path for plotting the edge in the base direction.  */
    pst_gr_mark_t emark;        /* Used for connected component enumeration. */
  } pst_gr_edge_data_t;
  /* Attributes of an (undirected) edge of the graph.
  
    The {path} defines a curve for plotting the edge. The curve starts
    at {org[0]}, goes near the vertices of {path}, and ends at {org[1]}.
    
    Let {ai} be the base arc of the edge. The field {enext[0]} is the
    arc that follows {ai} in the list of arcs that leave the origin
    vertex {org[0]}, in counterclockwise order. The field {enext[1]} is
    the arc that follows {pst_gr_arc_sym(ai)} in the list of arcs that leave
    the destination vertex {org[1]}, in counterclockwise order. The
    fields {eprev[0]} and {eprev[1]} are the same as {enext[0]} and
    {enext[1]}, except that the order is clockwise instead of
    counterclockwise.
    
    The {emark} is used internally when scanning the graph, and defines
    the vertex color when plotting. */

typedef struct pst_gr_t
  { uint32_t NV;                  /* Actual number of of vertices. */
    uint32_t NE;                  /* Actual number of (undirected) edges. */
    uint32_t NV_max;              /* Allocated space for vertices */
    uint32_t NE_max;              /* Allocated space for edges */
    pst_gr_vertex_data_t *vdata;  /* Vertex data, indexed {0..NV-1}*/
    pst_gr_edge_data_t *edata;    /* Edge data, indexed {0..NE-1} */
    /* Connection to height map: */
    uint32_t NX;                  /* Col count of height map. */
    uint32_t NY;                  /* Row count of height map. */
    pst_gr_vertex_t *vix;         /* Maps height map indices to vertex indices. */
  } pst_gr_t;
  /* The vertex data records are {vdata[0..NV-1]}, and the edge data records are
    {edata[0..NE-1]}. 
    
    If {NX} and {NY} are positive, the graph vertices are assumed to 
    correpond to a subset of the elements of an image (usually a height map) 
    with {NX} cols and {NY} rows.  In that case, the fields {x} and {y}
    of each vertex data record are the col and row indices of the correspondng 
    pixel, in the ranges {0..NX-1} and {0..NY-1}.  Also in that case,
    the table {vix} will have {NX*NY} elements, and {vix[NX*y + x]}
    will be the index of the vertex corresponding to pixel {[x,y]};
    or {pst_GR_NONE} if there is no such vertex. */

pst_gr_t *pst_gr_new(uint32_t NV_max, uint32_t NE_max, uint32_t NX, uint32_t NY);
  /* Creates a new graph record {gr}. Allocates the lists {gr.vdata},
    {gr.edata}, and {gr.hedge} with sizes {NV_max}, {NE_max}, and {NE_max}
    but initializes {gr.NV} and {gr.NE} to zero.  The {NX} and {NY}
    are both positive, allocates the table {gr.vix} with size {NX*NY},
    filling it with {pst_gr_NONE}. */

void pst_gr_free(pst_gr_t *gr);
  /* Reclaims all storage used by the graph {gr}, including the vertex
    and edge data records, the {gr.vix} table (if any) and the plot
    paths. */

pst_gr_edge_t pst_gr_arc_edge(pst_gr_arc_t ai);
  /* Returns the index of the undirected edge of the mesh 
    underlying the arc {ai}. If {ai} is {pst_gr_NONE}, returns {pst_gr_NONE}. */

pst_gr_dir_bit_t pst_gr_arc_dir_bit(pst_gr_arc_t ai);
  /* Returns the direction bit (0 or 1) of of the arc {ai} 
    (which must not be {pst_gr_NONE}). */
    
pst_gr_arc_t pst_gr_orient_edge(pst_gr_edge_t ei, pst_gr_dir_bit_t db);
  /* Returns the arc on the (unoriented) edge {ei} with direction bit {db}. */

pst_gr_edge_t pst_gr_arc_sym(pst_gr_arc_t ai);
  /* The arc with same edge s {ai} but opposite orientation. */

pst_gr_vertex_t pst_gr_arc_org(pst_gr_t *gr, pst_gr_arc_t ai);
  /* Returns the index of the vertex that is the origin of the arc {ai}
    (which must not be {pst_gr_NONE}.  Namely, {dr.edata[ei].org[db]}
    where {ei} is {pst_gr_arc_edge(ai)} and {db} is {pst_gr_arc_dir_bit(ai)}. */

pst_gr_vertex_t pst_gr_arc_dst(pst_gr_t *gr, pst_gr_arc_t ai);
  /* Returns the index of the vertex that is the destination of the arc {ai}
    (which must not be {pst_gr_NONE}.  Namely, {dr.edata[ei].org[1-db]}
    where {ei} is {pst_gr_arc_edge(ai)} and {db} is {pst_gr_arc_dir_bit(ai)}. */

pst_gr_arc_t pst_gr_arc_onext(pst_gr_t *gr, pst_gr_arc_t ai);
  /* The arc with same origin vertex {org} as {ai} that follows {ai}
    (which must not be {pst_gr_NONE}) in counterclockwise order around
    {org}. Namely, {dr.edata[ei].enext[db]} where {ei} is
    {pst_gr_arc_edge(ai)} and {db} is {pst_gr_arc_dir_bit(ai)}. */

pst_gr_arc_t pst_gr_arc_oprev(pst_gr_t *gr, pst_gr_arc_t ai);
  /* The arc with same origin vertex {org} as {ai} that precedes {ai}
    (which must not be {pst_gr_NONE}) in counterclockwise order around
    {org}. Namely, {dr.edata[ei].enext[db]} where {ei} is
    {pst_gr_arc_edge(ai)} and {db} is {pst_gr_arc_dir_bit(ai)}. */

pst_gr_arc_t pst_gr_arc_lnext(pst_gr_t *gr, pst_gr_arc_t ai);
  /* The arc with same left face {lft} as {ai} (which must not be
    {pst_gr_NONE}) that follows {ai} in counterclockwise order around
    {lft}.  namely {pst_gr_arc_oprev(gr, pst_gr_arc_sym(ai))}. */

double pst_gr_arc_weight(pst_gr_t *gr, pst_gr_arc_t ai);
  /* Gets the weight of the (undirected) edge underlying arc {ai} (which
    must not be {pst_gr_NONE}). The result is equal to
    {pst_gr_arc_weight(gr,pst_gr_arc_sym(ai))}. */

double pst_gr_arc_delta(pst_gr_t *gr, pst_gr_arc_t ai);
  /* Gets the supposed height difference (/delta/) along the (directed)
    arc {ai} (which must not be {pst_gr_NONE}). The result for
    {pst_gr_arc_sym(ai)} will be the negative of this. */

pst_gr_path_t pst_gr_arc_path(pst_gr_t *gr, pst_gr_arc_t ai);
  /* Gets the plotting path associated with the directed arc {ai} (which
    must not be {pst_gr_NONE}). */

r2_t pst_gr_arc_start_dir(pst_gr_t *gr, pst_gr_arc_t ai);
  /* Returns the unit vector {ui} that is the direction of the path of
    of the arc {ai} (which must not be {pst_gr_NONE}) as it leaves its
    origin vertex. */
    
uint32_t pst_gr_outdegree(pst_gr_t *gr, pst_gr_vertex_t vi);
  /* Returns the number of arcs out of the vertex with index {vi} (which
    must not be {pst_gr_NONE}). */

pst_gr_arc_t pst_gr_get_connecting_arc(pst_gr_t *gr, pst_gr_vertex_t vi0, pst_gr_vertex_t vi1);
  /* Returns the arc that goes from vertex number {vi0} to vertex number {vi1},
    of {pst_gr_NONE} if there is no such arc. */
  
pst_gr_vertex_t pst_gr_find_nearest_vertex(pst_gr_t *gr, r2_t *p);
  /* Returns the index of the vertex of {gr} that is closest to {p}.
    Fails if {gr} has no vertices. */

pst_gr_arc_t pst_gr_find_enclosing_face(pst_gr_t *gr, int32_t x, int32_t y);
  /* The graph {gr} must have an associated image (meaning that {gr->NX}
    and {gr->NY} must be positive). Returns and arc index {ai} such that
    the left face of {ai} contains the center of the map's grid cell
    with corners {(x,y)} and {(x+1,y+1)}.
    
    !!! Currently an incorrect hack !!! */

pst_gr_vertex_t pst_gr_get_vertex_from_map_indices(pst_gr_t *gr, int32_t x, int32_t y);
  /* The graph {gr} must have an associated map (meaning {gr->NX} and {gr->NY} must be positive).
    Returns a vertex {id} given the indices {x} and {y} of the corresponding pixel
    in that map.  The indices must be in {0..gr->NX-1}
    and {0..gr->NY-1}, respectively.  If there is no such vertex,
    returns {pst_gr_NONE}. */

void pst_gr_get_map_indices_from_vertex(pst_gr_t *gr, pst_gr_vertex_t vi, int32_t *x_P, int32_t *y_P);
  /* The graph {gr} must have an associated map (meaning {gr->NX} and {gr->NY} must be positive).
    If vertex {vi} (which must not be {NONE}) has has a corresponding pixel in
    that map, returns the indices in {*x_P} and {*y_P} the col and row indices of that pixel
    (which will be in {0..gr->NX-1} and {0..gr->NY-1}, respectively). 
    If vertex {vi} has no corresponding pixel, returns {*x = *y = -1}. */

pst_gr_vertex_t pst_gr_add_vertex
  ( pst_gr_t *gr,
    int32_t x,
    int32_t y,
    r2_t coords 
  );
  /* Adds to {gr} a new vertex with plot coordinates {coords},
    incrementing {gr.NV}. 
    
    If the vertex corresponds to a pixel in some associated image (such
    as a height map), the parameters {x,y} are the indices of that
    pixel, and must be in the ranges {0..gr->NX-1} and {0->gr->NY-1}. If
    there is no associated image, or the vertex is not associated with
    any pixel, both {x} and {y} should be {-1}.
    
    On exit, the new vertex will have index {vi = gr.NV-1}. The arc
    {gr.vdata[vi].aout} will be set to {pst_gr_NONE}. */

pst_gr_arc_t pst_gr_add_edge
  ( pst_gr_t *gr,
    pst_gr_vertex_t org,
    pst_gr_vertex_t dst,
    double d,
    double w,
    pst_gr_path_t P,
    bool_t setLinks
  );
  /* Adds to {gr} a new edge, incrementing {gr.NE}. On exit, the new
    edge will have index {ei = gr.NE-1}. The attributes {ed.org[0..1]}
    of the edge data record {ed=gr.edata[ei]} will be the vertices with
    indices {org} abd {dst}, respectively. The fields {ed.delta},
    {ed.weight, and {ed.path} will be set to {d,w,P}.
    
    Returns the base arc {ai} on the edge (oriented from {org} to {dst}).
    
    If {setLinks} is false, sets the links {ed.oprev[0..1]} and
    {ed.onext[0..1]} {pst_gr_NONE}, and does not change the {.oprev} and
    {.onext} of any other edges or the {.aout} field of any vertex. This
    option is to be used when those fields are provided through some
    other process.
    
    If {setLinks} is true, modifies the fields {eprev} and {enext} of
    the new edge and possibly other edges so that, on exit, the cycles
    defined by those links are consistent with the geometry of the
    edges. Namely, let {ai} be the base arc of the new edge {ei} upon
    exit. If on entry there are already arcs leaving vertex {org}, let
    {aj} and {ak} be arcs out of {org} such that the vector
    {pst_gr_arc_dir_vector(gr,ai)} lies between
    {pst_arc_dir_vector(gr,aj)} and {pst_arc_dir_vector(gr,ak)} in
    counterclockwise order around {org}. Then, on exit,
    {pst_gr_arc_onext(gr,aj)} and {pst_gr_arc_oprev(gr,ak)} will be
    {ai}. If on entry there are no arcs out of {org}, then on exit we
    will have {pst_gr_arc_onext(gr,ai)} and {pst_gr_arc_oprev(gr,ai)}
    both equal to {ai}; and the {.aout} field of {org} will be {ai}. The
    analogous condition will be true for the arc {pst_gr_arc_sym(ai)}
    and the vertex {dst}. This option requires that the {enext},
    {eprev}, and {aout} be consistent on entry, an dthey are left
    consistent on exit. */

void pst_gr_check_consistency(pst_gr_t *gr);
  /* Runs some validity checks on the graph {gr}.  Bombs out with messages
    if they fail. */
  
bool_t pst_gr_equal(pst_gr_t *gr, pst_gr_t *hr);
  /* Compares the two graphs {gr} and {hr}. If they are isomprphic and 
    isometric, returns {TRUE}.  If there are substantial
    differences, writes them to {stderr} and returns {FALSE}. */

double pst_gr_compute_left_face_curl(pst_gr_t *gr, pst_gr_arc_t ai);
  /* Estimates the net integral of the curl of the gradient field inside 
    the left face of arc {ai}, as the net sum of the deltas of the arcs
    in that face. */

r2_t pst_gr_left_face_barycenter(pst_gr_t *gr, pst_gr_arc_t ai);
  /* Computes the barycenter of the left face of arc {ai}, as the simple average of 
    the coordinates of its vertices.  Currently ignores the edge paths. */

void pst_gr_arc_print(FILE *wr, pst_gr_t *gr, pst_gr_arc_t ai);
  /* Writes to {wr} the arc {ai} of graph {gr}.  If {ai} is 
    {pst_gr_NONE}, prints "-1", else prints "{ei}:{db}" where 
    {ie} is the (undirected) edge index of {ai}
    and {db} is its direction bit. */
    
/* !!! There should be {pst_gr_MAX_VERTICES,pst_gr_MAX_EDGES} */

#define pst_gr_MAX_IMG_SIZE 4096
  /* Should be enough. */

#define pst_gr_FILE_TYPE "pst_gr_t"
  /* File type for {pst_gr_write_file} and {pst_gr_read_file}. */

#define pst_gr_FILE_VERSION "2025-03-12"
  /* File version for {pst_gr_write_file} and {pst_gr_read_file}. */

#endif


