/* Last edited on 2025-04-27 12:05:06 by stolfi */
/* Graph encoding of gradient maps for slope-to-height integration. */

#ifndef psgr_H
#define psgr_H

/* Improved version of {pst_graph.h} that emulates {haf.h} for the topology info. */ 
/* Created by Rafael F. V. Saracchini */

#include <stdint.h>

#include <r2.h>
#include <i2.h>
#include <psgr_types.h>
#include <psgr_path.h>

/* This interface defines a data structure to represent a planar graph properly embeded 
  in the plane, where each directed edge specifies the difference between
  the heights of the origina and destination vertices, and a reliability weight 
  for the same. */

psgr_t *psgr_new(uint32_t NV_max, uint32_t NE_max, uint32_t NX, uint32_t NY);
  /* Creates a new graph record {gr}. Allocates the lists of vertices
    and edges with initial sizes {NV_max} and {NE_max} but initializes
    the vertex and edge sets to empty. The underlying image will be
    assumed to have have {NX} cols and {NY} rows (which must be
    positive), so that vertex pixel indices will range in
    {{0..NX-1}×{0..NY-1}}. */

void psgr_free(psgr_t *gr);
  /* Reclaims all storage used by the graph {gr}, including internal tables and 
    edge paths. */
    
uint32_t psgr_num_vertices(psgr_t *gr);
  /* The number of vertices in {gr} (including those marked as deleted). */
    
uint32_t psgr_num_edges(psgr_t *gr);
  /* The number of edges in {gr} (including those marked as deleted). */

i2_t psgr_image_size(psgr_t *gr);
  /* Returns the number of cols and rows of the image associated with graph {gr},
    as  defined in {psgr_new}, packaged as an {i2_t} record. */

/* TOPOLOGY OPERATORS */

/* Unless said otherwise, every arc, edge, or vertex argument 
  may not be null ({psgr_arc_NONE}, {psgr_edge_NONE}, {psgr_vertex_NONE}). */

psgr_edge_t psgr_arc_edge(psgr_arc_t a);
  /* Returns the index of the undirected edge of the mesh 
    underlying the arc {a}. If {a} is {psgr_arc_NONE},
    returns {psgr_edge_NONE}. */

psgr_dir_bit_t psgr_arc_dir_bit(psgr_arc_t a);
  /* Returns the direction bit (0 or 1) of of the arc {a}. */

psgr_arc_t psgr_orient_edge(psgr_edge_t e, psgr_dir_bit_t db);
  /* Returns the arc on the (unoriented) edge {e} with direction bit {db}. */

void psgr_arc_split
  ( psgr_arc_t a,
    psgr_edge_t *e_P,
    psgr_dir_bit_t *db_P
  );
  /* Given an arc {a}, obtains its underlying undirected edge {e} and its
    direction bit {db}.  These values
    are returned in {*e_P} and {*db_P}, if they are not
    {NULL}.  Fails if {a} is {psgr_arc_NONE}. */

psgr_arc_t psgr_arc_sym(psgr_arc_t a);
  /* The arc with same edge s {a} but opposite orientation. */

psgr_vertex_t psgr_arc_org(psgr_t *gr, psgr_arc_t a);
  /* Returns the index of the vertex that is the origin of the arc {a}. */

psgr_vertex_t psgr_arc_dst(psgr_t *gr, psgr_arc_t a);
  /* Returns the index of the vertex that is the destination of the arc {a}. */

psgr_arc_t psgr_arc_onext(psgr_t *gr, psgr_arc_t a);
  /* The arc with same origin vertex {org} as {a} that follows {a}
    in counterclockwise order around {org}. */

psgr_arc_t psgr_arc_oprev(psgr_t *gr, psgr_arc_t a);
  /* The arc with same origin vertex {org} as {a} that precedes {a}
    in counterclockwise order around {org}. */

psgr_arc_t psgr_arc_dnext(psgr_t *gr, psgr_arc_t a);
  /* The arc with same destination vertex {dst} as {a} that follows {a}
    in counterclockwise order around {dst}.
    Namely {psgr_arc_sym(psgr_arc_onext(gr,psgr_arc_sym(a)))}. */

psgr_arc_t psgr_arc_dprev(psgr_t *gr, psgr_arc_t a);
  /* The arc with same destination vertex {dst} as {a} that precedes {a}
    in counterclockwise order around {dst}.
    Namely {psgr_arc_sym(psgr_arc_oprev(gr,psgr_arc_sym(a)))}. */

psgr_arc_t psgr_arc_lnext(psgr_t *gr, psgr_arc_t a);
  /* The arc with same left face {lft} as {a} that follows {a} 
    in counterclockwise order around {lft}.  
    Namely {psgr_arc_oprev(gr, psgr_arc_sym(a))}. */

psgr_arc_t psgr_arc_lprev(psgr_t *gr, psgr_arc_t a);
  /* The arc with same left face {lft} as {a} that precedes {a}
    in counterclockwise order around {lft}.  
    Namely {psgr_arc_sym(psgr_arc_onext(gr, a))}. */
    
psgr_arc_t psgr_vertex_out_arc(psgr_t *gr, psgr_vertex_t v);
  /* Returns one arc of the onext ring of {v}, or {psgr_arc_NONE}
    if {v} has no incident edges. */

uint32_t psgr_vertex_outdegree(psgr_t *gr, psgr_vertex_t v);
  /* Returns the number of arcs out of the vertex with index {v}. */

psgr_arc_t psgr_get_connecting_arc(psgr_t *gr, psgr_vertex_t vi0, psgr_vertex_t vi1);
  /* Returns the arc that goes from vertex number {vi0} to vertex number {vi1},
    of {psgr_arc_NONE} if there is no such arc. */

void psgr_arc_print(FILE *wr, psgr_t *gr, psgr_arc_t a);
  /* Writes to {wr} the arc {a} of graph {gr}.  If {a} is 
    {psgr_arc_NONE}, prints "-1", else prints "{e.ix}:{db}" where 
    {e} is the (undirected) edge of {a} and {db} is its direction bit. */
    
/* ARC/EDGE/VERTEX ATTRIBUTES */

double psgr_arc_delta(psgr_t *gr, psgr_arc_t a);
  /* Gets the supposed height difference (/delta/) along the (directed)
    arc {a} (which must not be {psgr_arc_NONE}). The result for
    {psgr_arc_sym(a)} will be the negative of this. */

void psgr_arc_set_delta(psgr_t *gr, psgr_arc_t a, double d);
  /* Sets the delta of arc {a} to {d}.  Equivalent to 
    {psgr_arc_set_delta(gr,psgr_arc_sym(a), -d)} (note the sign). */

double psgr_arc_weight(psgr_t *gr, psgr_arc_t a);
  /* Gets the weight of the (undirected) edge underlying arc {a} (which
    must not be {psgr_arc_NONE}). The result is equal to
    {psgr_arc_weight(gr,psgr_arc_sym(a))}. */

void psgr_arc_set_weight(psgr_t *gr, psgr_arc_t a, double w);
  /* Sets the weight of arc {a} to {w}.  Equivalent to 
    {psgr_arc_set_weight(gr,psgr_arc_sym(a),w)}. */

psgr_mark_t psgr_vertex_mark(psgr_t *gr, psgr_vertex_t v);
psgr_mark_t psgr_edge_mark(psgr_t *gr, psgr_edge_t e);
  /* Get the current mark of the vertex {v} and edge {e}, respectively. */

void psgr_vertex_set_mark(psgr_t *gr, psgr_vertex_t v, psgr_mark_t mrk);
void psgr_edge_set_mark(psgr_t *gr, psgr_edge_t e, psgr_mark_t mrk);
  /* Set the current mark of the vertex {v} and edge {e}, respectively, to {mrk}. */

/* GRAPH GEOMETRY */

i2_t psgr_vertex_pixel(psgr_t *gr, psgr_vertex_t v);
  /* If vertex {v} (which must not be {psgr_vertex_NONE}) has has a
    corresponding pixel in the underlying image, returns the (integer)
    col and row indices of that pixel. They will be in {0..NX-1} and
    {0..NY-1}, respectively, where {NX} and {NY} are the col and row
    counts of the image. If vertex {v} has no corresponding pixel,
    returns {*x = *y = -1}. */
    
psgr_path_t psgr_arc_path(psgr_t *gr, psgr_arc_t a);
  /* Gets the plotting path associated with the directed arc {a} (which
    must not be {psgr_arc_NONE}). */

r2_t psgr_arc_starting_dir(psgr_t *gr, psgr_arc_t a);
  /* Returns the unit vector {ui} that is the direction of the path of
    of the arc {a} (which must not be {psgr_arc_NONE}) as it leaves its
    origin vertex. */

double psgr_arc_left_face_curl(psgr_t *gr, psgr_arc_t a);
  /* Estimates the net integral of the curl of the gradient field inside 
    the left face of arc {a}, as the net sum of the deltas of the arcs
    in that face. */

r2_t psgr_arc_left_face_barycenter(psgr_t *gr, psgr_arc_t a);
  /* Computes the barycenter of the left face of arc {a}, as the simple average of 
    the pixel indices of its vertices.  Currently ignores the edge paths. */
  
psgr_vertex_t psgr_find_nearest_vertex(psgr_t *gr, r2_t *p);
  /* Returns the index of the vertex of {gr} that is closest to {p}.
    Fails if {gr} has no vertices. */

psgr_arc_t psgr_find_enclosing_face(psgr_t *gr, r2_t *p);
  /* The graph {gr} must have an associated image (meaning that {gr->NX}
    and {gr->NY} must be positive). Returns and arc index {a} such that
    the left face of {a} contains the point {p}.  If there is no such 
    face, returns {psgr_arc_NONE}.
    
    !!! Currently an incorrect hack !!! */

/* VERTEX-PIXEL MAPPING */

psgr_vertex_t psgr_vertex_from_pixel(psgr_t *gr, i2_t *ip);
  /* Returns a vertex {v} given its pixel indices {ip.c[0],ip.c[1]}.  The indices must be in {0..NX-1}
    and {0..NY-1}, respectively, where {NX} and {NY} are the col and row counts of the 
    underlying map.  If there is no such vertex,  returns {psgr_vertex_NONE}. */

/* GRAPH CREATION AND EXPANSION */

psgr_vertex_t psgr_add_vertex(psgr_t *gr, i2_t *ip);
  /* Adds to {gr} a new vertex {v} with pixel indices {[ip.c[0],ip.c[1]]},
    which must be in the ranges {0..NX-1} and {0..NY-1},
    where {(NX,NY)} is {psgr_image_size(gr)}.  
    
    The vertex index {v.ix} will be {NV-1} where {NV} is
    {psgr_num_vertices(gr)} on exit. The vertex mark will be set to
    {psgr_mark_CLEAR} and it will have no incident edges. */

psgr_arc_t psgr_add_edge
  ( psgr_t *gr,
    psgr_vertex_t org,
    psgr_vertex_t dst,
    double d,
    double w,
    psgr_path_t P
  );
  /* Adds to {gr} a new edge {e} connecting vertices {org} and {dst}.
    Returns the arc {a} of that edge that is oriented from {org} to
    {dst}.
  
    The new edge index {e.ix} will be {NE-1} where {NE} is
    {psgr_num_edges(gr)} on exit. The delta, weight, and path of the
    arc {a} will be set to {d,w,P}.
    
    The procedure also inserts {a} and its symmetrical into the
    onext rings of {org} and {dst}, respectively, based on the initial
    directions of the respective oriented paths.  This will define the
    values of {psgr_arc_onext(a)}, {psgr_arc_oprev(a)}, 
    {psgr_arc_dnext(a)}, and {psgr_arc_dprev(a)} to be arcs of those
    two rings. */

/* VERTEX AND EDGE REMOVAL */

void psgr_edge_remove(psgr_t *gr, psgr_edge_t e, int32_t indent);
  /* The edge {e} must not be {psgr_edge_NONE} and the {.emark} of
    its data record must be {psgr_mark_CLEAR}. Disconnects the edge
    {e} from the graph, by removing both its arcs from the respective
    {onext} rings and vertex {.aut} pointers, then making those arcs
    into singleton rings. The {.emark} field is then set to
    {psgr_mark_DELETED}. The {.path} record, if any, is not 
    reclaimed. The {.weight} and {.data} fields are not changed.
    
    It {indent} is non-negative, prints diagnostics indented by {indent} spaces. */

void psgr_vertex_remove(psgr_t *gr, psgr_vertex_t v, int32_t indent);
  /* Effectively removes the vertex {v} from the graph {gr}
    by setting its {.vmark} field to {psgr_mark_DELETED}.
    The (which must not be from the graph {gr}.
    The vertex {v} must not be {psgr_vertex_NONE}
    and must be isolated, i.e. it must not be origin or destination
    of any edge (namely its {.aout} must be {psgr_arc_NONE}).
    Its {vmark} cannot be {psgr_mark_DELETED} or {psgr_mark_KEEP}.
    
    It {indent} is non-negative, prints diagnostics indented by {indent} spaces. */ 
    
/* DEBUGGING TOOLS */

void psgr_check_consistency(psgr_t *gr);
  /* Runs some validity checks on the graph {gr}.  Bombs out with messages
    if they fail.
    
    Every vertex mark must be in {{CLEAR,KEEP,KILL,DELETE}}, and every
    edge mark must be in {{CLEAR,DELETED}}. Morever, edge incident to a
    vertex marked {DELETED} mus be marked {DELETED} too. */
  
bool_t psgr_equal(psgr_t *gr, psgr_t *hr);
  /* Compares the two graphs {gr} and {hr}. If they are isomprphic and 
    isometric, returns {TRUE}.  If there are substantial
    differences, writes them to {stderr} and returns {FALSE}. */

#endif


