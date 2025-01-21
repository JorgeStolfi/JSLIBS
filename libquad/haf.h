#ifndef haf_H
#define haf_H

/* The half-edge data structure to encode the topology of 2D meshes. */
/* Last edited on 2025-01-10 08:14:31 by stolfi */

#include <stdint.h>

#define haf_H_copyright \
  "Copyright © 2023 State University of Campinas (UNICAMP).\n\n" jslibs_copyright

#include <jslibs_copyright.h>
#include <bool.h>
#include <vec.h>

/* THE HALF-EDGE STRUCTURE 

  The half-edge data structure is used to encode the topology of a
  two-dimensional /mesh/, namely a subdivision of a two-dimensional
  compact topological manifold (the mesh's /surface/) into a finite
  number of /elements/, which are open topological disks (the /faces/)
  open topological segments (the /edges/) and points (the /vertices/),
  with certain incidence and adjacency conditions. 
  
  The manifold is supposed to be oriented, meaning that at each point
  there is a defined "positive" or "counterclockwise" sense of turning,
  that can be continuously transported over each connected component of
  the manifold without contradiction.
  
  An /arc/ of the mesh is an edge in a specific orientation. The two
  arcs of the same edge are said to be /opposite/ of each other. The
  (/oriented/) /border/ of each face is a circular list of alternating
  arcs and vertices, where each arc is oriented in the direction that
  matches counterclockwise turning on the face. Each arc may appear only
  once in this list, and only on one face (but an arc and its opposite may
  both appear on the same face).
  
  The half-edge structure consists of /arc records/, each representing,
  by definition, one arc of the mesh. The structure defines two
  fundamental /walking operators/ on ac records, {.sym} and {.lnext}.
  The {.sym} operator is an involution of set of all arc records, with
  no fixed points: {a.sym != a} and {a.sym.sym = a} for every arc {a}. The
  records {a} and {a.sym} are interpreted as the oppositely oriented
  arcs on the same edge of the mesh.
  
  The {.lnext} operator is an arbirary permutation of the arc records,
  meaning that {a.lnext != b.lnext} if {a != b}. Each cycle of this
  permutation is interpreted as the list of arcs on the border of a face
  of the mesh.
  
  The {.oprev} operator is defined as the composition {.sym.lnext},
  applied in that order. Each cycle of {.oprev} corresponds to a unique
  vertex of the mesh, and consists of the arcs that leave that vertex,
  in /clockwise/ order around the vertex. Again, an arc {a} may appear
  only once in only one of these cycles; but {a} and {a.sym} may both
  appear on the same cycle, signifying that the corresponding edge is a
  loop with both endpoints on the same vertex.
  
  Several other walking operators can be defined in terms of these, such as
  {.lprev}, {.onext}, {.rnext}, {.rprev}, {.dprev}, {.dnext}. See below.
  
  Each arc record {a} also contains a non-negative integer /arc id/
  {a.aid} that can be chosen by the client of this interface. The arc
  ids of {a} and {a.sym} are always {2*a.eid} and {2*a.eid+1} or
  vice-versa, where {a.eid=a.sym.eid} is a client-chosen id number for
  the underlying undirected edge. The integers {a.aid} and {a.eid} can
  be used as indices in external tables to attach arbitrary properties
  to each arc or edhe of the mesh.
  
  The half-edge data structure proper has no explicit records for faces
  and vertices. Clients of this interface who need to attach properties
  to those elements of the mesh (such as vertex coordinates or face
  geometry) may define their own data records for them. In that case,
  they should probably set up tables {left[a.aid]} that gives the record
  for the face that has {a} on its border, and/or {org[a.aid]} for the
  record of the vertex that is the origin of arc {a}. See {haf_enum_faces}
  and {haf_enum_vertices} below. */

typedef struct haf_rep_t *haf_arc_t;
  /* An /arc reference/, a pointer to a record of the half-edge 
    structure representing an /arc/ (directed edge) of the mesh. */

vec_typedef(haf_arc_vec_t, haf_arc_vec, haf_arc_t);
  /* An extensible vector of arc references. */

/* FAST WALKING OPERATORS */

haf_arc_t haf_sym(haf_arc_t a);
  /* Same edge in opposite orientation.  Takes constant time. */
  
haf_arc_t haf_lnext(haf_arc_t a);
  /* Next arc /counterclockwise/ with the same left face. Takes constant time. */
 
haf_arc_t haf_rnext(haf_arc_t a);
  /* The next arc /counterclockwise/ with the same right face,
    namely {a.sym.lnext.sym}. Takes constant time. */
       
haf_arc_t haf_oprev(haf_arc_t a);
  /* The next arc /clockwise/ with the same origin vertex,
    namely {a.sym.lnext}. Takes constant time. */
  
haf_arc_t haf_dprev(haf_arc_t a);
  /* The next arc /clockwise/ with the same destination vertex
    namely {a.lnext.sym}. Takes constant time. */

/* SLOW WALKING OPERATORS */

haf_arc_t haf_lprev(haf_arc_t a);
  /* Inverse of {.lnext}: the next arc /clockwise/ with the same left face.
    Takes time proportional to the degree of that face. */

haf_arc_t haf_rprev(haf_arc_t a);
  /* Inverse of {.rnext}: next arc /clockwise/ with same right face.
    Takes time proportional to the degree of that face. */
      
haf_arc_t haf_onext(haf_arc_t a);
  /* Inverse of {.oprev}: the next arc /counterclockwise/ with the same origin vertex.
    Takes time proportional to the degree of that vertex. */
   
haf_arc_t haf_dnext(haf_arc_t a);
  /* Inverse of {.dprev}: the next arc /counterclockwise/ with the same destination vertex.
    Takes time proportional to the degree of that vertex. */

/* ARC DIRECTION BIT 

  The /direction bit/ of an arc, either 0 or 1, distinguishes its from
  its opposite. Thus {a} and {a.sym} have direction bits 0 and 1 or
  vice-versa. An arc with direction bit 0 is said to be the /base arc/
  of the underlying edge. */

typedef uint8_t haf_dir_bit_t;
   /* A data type used to hold a few bits (usually at the low-order end). */

haf_dir_bit_t haf_dir_bit(haf_arc_t a);
  /* The direction bit of arc {a}. */

/* UNORIENTED EDGES */

typedef struct haf_edge_rec_t *haf_edge_t;
  /* Reference to an undirected and unoriented edge {e}
    of the mesh.  It is shared by the two arcs that consist of {e}
    taken with both directions. */

vec_typedef(haf_edge_vec_t, haf_edge_vec, haf_edge_t);
  /* An externsible vector od {half_edge_t} references. */

haf_edge_t haf_edge(haf_arc_t a);
  /* Obtains the edge reference of an arc reference {a}. Satisfies
    {a.sym.edge == a.edge}.  */

haf_arc_t haf_orient(haf_edge_t e, haf_dir_bit_t db);
  /* Returns the arc {a} that has {haf_edge(a) = e} and {had_dir_bit(a) = db}.  */

haf_arc_t haf_base_arc(haf_edge_t e);
  /* The base arc of the edge {e} (the one with direction bit 0). */

/* EDGE IDENTIFIERS 
  
  The structure stores a identifying non-negative /id number/ for each
  undirected edge of the mesh. See {haf_set_edge_id} and
  {haf_renumber_edges} below for how that number gets defined.
  
  The edge identifier is often meant to identify each edge uniquely,
  e.g. for debugging or table indexing; but this interface does not 
  guarantee that by itself. */

typedef uint64_t haf_edge_id_t;  /* An edge ID number. */

#define haf_edge_id_MAX (haf_arc_id_MAX/2)
  /* Max edge identifier (so that the arc identifier does not overflow). */

haf_edge_id_t haf_edge_id(haf_edge_t e);
  /* Returns the identifier of the undirected edge {e},
    namely {aid/2} where {aid} is the id of the edge's base arc. */

typedef uint64_t haf_edge_count_t;  /* A count of edges in a data structure, table, etc.. */

#define haf_edge_count_MAX (haf_arc_count_MAX/2)
  /* Max number of edges in a structure. */

vec_typedef(haf_edge_id_vec_t, haf_edge_id_vec, haf_edge_id_t);
  /* An extensible vector of edge ids. */

/* ARC IDENTIFIERS 

  The /arc identifier/ of an arc {a} is {2*eid + db} where
  {eid} is the identifying number of the undirected edge, 
  and {db} is the direction bit. */

typedef uint64_t haf_arc_id_t;   /* An arc ID number. */

haf_arc_id_t haf_arc_id(haf_arc_t a);
  /* Returns the identifier of arc {a}. */

#define haf_arc_id_MAX (UINT64_MAX)
  /* Max arc identifier. */

typedef uint64_t haf_arc_count_t;   /* A count of arcs in a data structure, table, etc.. */

#define haf_arc_count_MAX (((uint64_t)1024)*1024*1024)
  /* Max number of arcs in a structure for sane table allocation. */

vec_typedef(haf_arc_id_vec_t, haf_arc_id_vec, haf_arc_id_t);
  /* An extensible vector of arc ids. */

/* BUILDING AND MODIFICATION */
  
haf_arc_t haf_make_stick(haf_edge_id_t eid);
  /* Creates a new half-edge structure consisting of a pair of
    {haf_arc_t} records, {a,b}, representing a mesh with sperical
    topology, one non-loop edge, one face, and two distinct vertices.
    The records will have {a.sym=b}, {b.sym=a}, {a.lnext=b},
    {b.lnext=a}. The arc ids will be {2*eid} and {2*eid+1}. */

haf_arc_t haf_make_loop(haf_edge_id_t eid);
  /* Creates a new half-edge structure consisting of a pair of
    {haf_arc_t} records, {a,b}, representing a mesh with spherical
    topology, one loop edge, one vertex, and two distinct faces. The
    records will have {a.sym=b}, {b.sym=a}, {a.lnext=a}, {b.lnext=b}.
    The arc ids will be {2*eid} and {2*eid+1}. */

void haf_edge_free(haf_edge_t e);
  /* Reclaims the space used by the edge {e}. 
  
    The edge must be isolated, meaning that {a.lnext} and {b.lnext} must
    be respectively {a} and {b} or {b} and {a}, where {a,b} are the two
    arcs on {e}. That is, {e} must be an isolated stick or an isolated
    loop.  The two arcs become invalid after this operation and should 
    not be used in any way. */

void haf_splice(haf_arc_t a, haf_arc_t b);
  /* Performs a splice (split and join) operation on the {.lnext} loops
    of the arcs {a} and {b} and their opposites. Namely, swaps {a.lnext}
    with {b.lnext}. 
    
    As a consequence, (1) if the left faces (the
    {.lnext} loops) of {a} and {b} were distinct, they becme one, and
    vice-versa; and (2) if the destination vertics ({.dprev} loops) of
    {a} and {b} were distinct, they become one. 
    
    This operation maintains the topological consistency of the data
    structure. It may however split a connected component of the surface
    into two separate components, or vice-versa; or may create or remove
    a handle on a component of the surface. In any case, the surface
    continues to be orientable, and the implied orientation of the new
    {a.left} and {b.left} face(s) is consistent with those of {a.right}
    and {b.right}, and that of any other face other than the original
    {a.left} and {b.left} is not affected. */
 
void haf_set_lnext(haf_arc_t a, haf_arc_t b);
  /* Sets {a.lnext} to be {b}.  This operation may break the topological consistency.
    The caller must make sure that eventually the consistency is restored, namely that 
    the {.lnext} operator is a permutation of all the {haf_arc_t} records. */

void haf_set_edge_id(haf_edge_t e, haf_edge_id_t eid);
  /* Sets the identifier of edge {e} to {eid}.  Thus, if {a} is the base arc of {e},
    it will set {haf_arc_id(a)} to {2*eid} and {haf_arc_id(a.sym)} to {2*eid+1}. */
  
/* DEBUGGING */

void haf_check_topology(haf_edge_count_t NE, haf_edge_t e[], haf_edge_id_t eid0, bool_t verbose);
  /* Expects that the vector {e[0..NE-1]} has the edges of the mesh, in
    order of edge id. Checks that the id of {e[ke]} is {eid0+ke} for
    all {ke} in {0..NE-1}. Checks that the {.lnext} operator is a
    permutation of the arcs derived from those edges, and that the
    {.sym} operator is an involution on those arcs, with no fixed
    point. */

#endif
