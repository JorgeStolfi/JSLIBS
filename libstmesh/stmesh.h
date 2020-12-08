/* Arbitrary triangle mesh, with quasimanifold topology. */
/* Last edited on 2016-05-03 15:45:18 by stolfilocal */

#ifndef stmesh_H
#define stmesh_H

#define _GNU_SOURCE
#include <stdio.h>

#include <vec.h>
#include <bool.h>
#include <i3.h>
#include <r2.h>
#include <r3.h>

/* TRIANGLE MESHES WITH TOPOLOGICAL INFORMATION */

typedef struct stmesh_rep_t *stmesh_t;
  /* An {stmesh_t} is a pointer to a data structure that represents a
    read-only list of triangles in {\RR^3}, with topological
    information. It has a list of distinct vertices, a list of distinct
    edges, and a list of distinct faces (triangles). The coordinates of
    each vertex are rounded to even multiples of a fundamental length
    unit {eps}, stored in the structure. */
     
float stmesh_get_eps(stmesh_t mesh);
  /* Returns the fundamental length unit (in mm) of the given mesh. */
    
void stmesh_get_bounding_box(stmesh_t mesh, i3_t *minQP, i3_t *maxQP);
  /* Stores into {*minQP,*maxQP} the lower and upper
    corner of the bounding box of {mesh}. */

void stmesh_print_bounding_box(FILE *wr, stmesh_t mesh);
  /* Writes to {wr} a readable display of the mesh's bounding box,
    in quantized coordinates and in {mm}. */

uint32_t stmesh_vert_count(stmesh_t mesh);
uint32_t stmesh_edge_count(stmesh_t mesh);
uint32_t stmesh_face_count(stmesh_t mesh);
  /* These procedures return the number of distinct vertices,
    edges, and faces in the given {mesh}, respectively. */

/* VERTICES, EDGES, AND FACES. */

/* The parts of the mesh - vertices, edges, and faces -- are represented
  by records with private format. The edge and face pointer types below
  may have their low-order bits modified to indicate an orientation for
  the element. Thus they must not be dereferenced directly. */

typedef struct stmesh_vert_rep_t *stmesh_vert_t;
  /* A pointer to a record that represents a vertex in a {stmesh_t}. */

typedef struct stmesh_edge_rep_t *stmesh_edge_t;
  /* A pointer to a record that represents an ORIENTED edge in a {stmesh_t};
    that is, an edge with a specific sense of traversal.  The orientation defines
    which endpoint is the /origin/ and which is the /destination/,
    and also what is the circular ordering of the faces around the edge, if there are
    more than two (by the right-hand rule). */
     
typedef struct stmesh_face_rep_t *stmesh_face_t;
  /* A pointer to a record that represents an ORIENTED triangle in a
    {stmesh_t}. The orientation specifies one of the two possible
    senses of traversal of its perimeter, and also designates a specific
    side as being the triangle's /base edge/. Thus, there are six distincg
    oriented triangles for each actual (unoriented) triangle. */

/* VERTEX OPERATIONS */    
 
i3_t stmesh_vert_get_pos(stmesh_vert_t v);
  /* Returns the quantized coordinates of the vertex {v}. */

r3_t stmesh_unround_point(i3_t *p, float eps);
  /* Converts the quantized point {*p} to a non-quantized point,
    with coordinates in mm, by multiplying it by {eps}. */

/* EDGE OPERATIONS */

stmesh_edge_t stmesh_edge_reverse(stmesh_edge_t e, int k);
  /* If {k} is odd, returns the oriented edge which is the same unoriented edge as
    {e}, but with the opposite orientation.  If {k} is even, retruns {e}
    itself. */

stmesh_edge_t stmesh_edge_natural(stmesh_edge_t e);
  /* Returns the same edge as {e}, in its "natural" orientation.
    That is, {stmesh_edge_natural(e)==stmesh_edge_natural(e')=e*}
    for any oriented edge {e}, where {e'} is the reverse of {e},
    and {e*} is either {e} or {e'}. */

stmesh_vert_t stmesh_edge_get_endpoint(stmesh_edge_t e, int k);
  /* Returns the vertex that is the origin ({k = 0}) or
    destination ({k = 1}) of {e}, taking its orientation into account. */

void stmesh_edge_get_endpoints(stmesh_edge_t e, stmesh_vert_t v[]);
  /* Stores into {v[0..1]} the end vertices of {e}, in the order determined
    by its orientation.  Namely {v[0]}  will be the origin, and {v[1]} the destintion. */
    
uint32_t stmesh_edge_degree(stmesh_edge_t e);
  /* Returns the degree of edge {e}, that is, the number of 
    faces incident to {e}. Ignores the orientation. */

/* FACE OPERATIONS */    

stmesh_edge_t stmesh_face_get_base(stmesh_face_t f);
  /* Return the oriented base edge of the oriented face {f}. */
  
stmesh_face_t stmesh_face_flip(stmesh_face_t f, int k);
  /* If {k} is odd, returns the oriented face which is the same
    unoriented face as {f}, but taken with the opposite orientation. If
    {k} is even, returns {f} itself.
    
    When the orientation is flipped, the base edge of the result will be
    the base edge of {f}, reversed. The new orientation specifies that
    the sides of {f} will be traversed in the opposite order, starting from
    that new base edge, and each side will traversed in the opposite
    sense. */

stmesh_face_t stmesh_face_shift(stmesh_face_t f, int k);
  /* Returns the oriented face that is the same unoriented
    face as {f}, with the perimeter traversed in the same sense, but whose
    base edge is {k} sides ahead of {f}'s base edge, in that sense of traversal.
    
    In particular, {stmesh_face_shift(f,k)} is the same as
    {stmesh_face_shift(f,k+3)}, and
    {stmesh_face_shift(stmesh_face_flip(f),k)} is
    {stmesh_face_flip(stmesh_face_shift(f,-k))}, for all {k}. */

stmesh_face_t stmesh_face_natural(stmesh_face_t f);
  /* Returns the oriented face that is the same unoriented
    face as {f}, in its "natural" sense and base edge.
    That is, {stmesh_face_natural(f)==stmesh_face_natural(f')=f*}
    for any oriented versions {f} and {f'} of the same 
    unoriented face; where {f*} is a particular
    flipped and/or shifted version of {f}. */

void stmesh_face_get_sides(stmesh_face_t f, stmesh_edge_t e[]);
  /* Returns oriented edges that are the sides of the 
    face {f}, starting with its base edge and progressing along the perimeter in the
    direction specified by its orientation.  The edges will be oriented so as to 
    agree with that sense of traversal. */

void stmesh_face_get_corners(stmesh_face_t f, stmesh_vert_t v[]);
  /* Stores in {v[0..2]} the three corners of the oriented face {f},
    starting with the vertex opposte to its base edge and progressing along the perimeter in the
    direction specified by its orientation. */

void stmesh_face_get_zrange(stmesh_face_t f, int32_t *minZP, int32_t *maxZP);
  /* Returns in {*minZP,*maxZP} the minimum and maximum quantized {Z} coordinates of the 
    corners of the triangle {f}.  Ignores the triangle's orientation. */
 
/* GEOMETRIC SLICING OPERATIONS */

bool_t stmesh_edge_crosses_plane(stmesh_edge_t e, int32_t pZ);
  /* Returns {TRUE} iff the edge {e} crosses the horizontal plane 
    with quantized {Z}-coordinate {pZ}. The outcome is undefined 
    if the plane contains one or both endpoints of the edge. */
 
r2_t stmesh_edge_plane_intersection(stmesh_edge_t e, int32_t pZ, double eps);
  /* Computes the coordinates (in mm) of the point where the edge {e}
    intersects the slicing plane with quantized {Z}-coordinate {pZ}.
    The outcome is undefined if the plane does not intersect the edge,
    or contains one or both endpoints of the edge. */

void stmesh_face_get_sliced_sides(stmesh_t mesh, stmesh_face_t f, int32_t pZ, stmesh_edge_t e[]);
  /* Given the quantized {Z}-coordinate {pZ} of a slicing plane, finds
    the two sides of the oriented face {f} that are crossed by the
    plane. Assumes that the plane intersects the face but does not go
    through any of its vertices. Returns the two edges in arbitrary order,
    in their natural orientation. */

/* ELEMENT INDEXING */

typedef uint32_t stmesh_vert_unx_t;
  /* An integer that identifies a
    vertex of a mesh, namely the index of its vertex record in 
    in the table {mesh.v}.  Must be in {0..nv-1}, where 
    {nv=stmesh_vert_count(mesh)}. */

typedef uint32_t stmesh_edge_unx_t;
  /* An integer that identifies an UNORIENTED edge in a
    {stmesh_t}, namely the index of its vertex record in in the table
    {mesh.e}. Must be in the range {0..ne-1}, where
    {ne=stmesh_edge_count(mesh)}. */

typedef uint32_t stmesh_face_unx_t;
  /* An integer that identifies an UNORIENTED triangle in a
    {stmesh_t}, namely the index of its vertex record in in the table
    {mesh.f}. Must be in the range {0..nf-1}, where {ne =
    stmesh_edge_count(mesh)}. */

stmesh_vert_t stmesh_get_vert(stmesh_t mesh, stmesh_vert_unx_t uxv);
stmesh_edge_t stmesh_get_edge(stmesh_t mesh, stmesh_edge_unx_t uxe);
stmesh_face_t stmesh_get_face(stmesh_t mesh, stmesh_face_unx_t uxf);
  /* Get the vertex, edge, or face of {mesh} with the speficied index.
    The index must be in the appropriate range. Edges and faces are
    returned in their natural orientation. */

stmesh_vert_unx_t stmesh_vert_get_unx(stmesh_t mesh, stmesh_vert_t v);
stmesh_edge_unx_t stmesh_edge_get_unx(stmesh_t mesh, stmesh_edge_t e);
stmesh_face_unx_t stmesh_face_get_unx(stmesh_t mesh, stmesh_face_t f);
  /* Returns the index of the vertex, edge, or face in {mesh}.
    disregards the orientation of edges and faces. */

/* CONSTRUCTION */

stmesh_t stmesh_new_desc(float eps, uint32_t nv_max, uint32_t ne_max, uint32_t nf_max);
  /* Creates a new descriptor record for a mesh, with capacity for
    {nv_max} vertices, {ne_max} edges, and {nf_max} faces,
    but with no edges, no vertices, and no faces. */

stmesh_vert_unx_t stmesh_add_vert(stmesh_t mesh, i3_t *pos);
  /* Adds a vertex to {mesh}, given its quantized coordinates {*pos}.
    Returns its index in the mesh. Updates the bounding box of the mesh
    to include that vertex. The vertex had better be used later in some
    edge. */

stmesh_edge_unx_t stmesh_add_edge(stmesh_t mesh, stmesh_vert_unx_t uxv[]);
  /* Adds an edge to the mesh, given the indices {uxv[0..1]} of its
    endpoint vertices. The indices must be in increasing order. Sets its
    degree to zero. Returns the index of the edge in the mesh. The edge
    had better be used later by some face. Assumes that {mesh->v} is set. */
        
stmesh_face_unx_t stmesh_add_face(stmesh_t mesh, stmesh_edge_unx_t uxe[]);
  /* Adds a face to the mesh, given the indices {uxe[0..2]} of the
    unoriented edges that bound the face. The indices should be sorted
    in increasing order. Increments the degree of the edges that are
    sides of {f}. Returns the index of the face in the mesh. Assumes
    that {mesh.{nv,v,ne,e}} are set. */

void stmesh_free(stmesh_t mesh);
  /* Releases all heap storage allocated for {mesh}, including the 
    descriptor record. */

/* MISCELLANEOUS */

void stmesh_print_edge_degrees(FILE *wr, stmesh_t mesh);
  /* Counts how many edges have degree {0,1,...}, 
    prints them to {wr}.  Aborts with error message 
    if there is an edge with degree 0. */ 

void stmesh_print(FILE *wr, stmesh_t mesh);
void stmesh_vert_print(FILE *wr, stmesh_t mesh, stmesh_vert_t v);
void stmesh_edge_print(FILE *wr, stmesh_t mesh, stmesh_edge_t e);
void stmesh_face_print(FILE *wr, stmesh_t mesh, stmesh_face_t f);
  /* Writes legibly to {wr} the address, number, and fields of the element. */

void stmesh_check(stmesh_t mesh);
void stmesh_vert_check(stmesh_t mesh, stmesh_vert_t v);
void stmesh_edge_check(stmesh_t mesh, stmesh_edge_t e);
void stmesh_face_check(stmesh_t mesh, stmesh_face_t f);
  /* Runs some consistency checks on the topological operations
    on those elements. */

/* MESH BUILDING FROM TABLES */

typedef struct stmesh_vert_unx_pair_t { stmesh_vert_unx_t c[2]; } stmesh_vert_unx_pair_t;
  /* A pair of vertex indices. */

typedef struct stmesh_edge_unx_triple_t { stmesh_edge_unx_t c[3]; } stmesh_edge_unx_triple_t;
  /* A triplet of unoriented edge indices. */

stmesh_t stmesh_build
  ( float eps, 
    uint32_t nv, 
    i3_t vpos[], 
    uint32_t ne, 
    stmesh_vert_unx_pair_t endv[],
    uint32_t nf, 
    stmesh_edge_unx_triple_t side[],
    bool_t checkSorted
  );
  /* Builds a mesh data structure with fundamental length unit {eps},
    {nv} vertices, {ne} unoriented edges, and {nf} unoriented faces,
    given the vertex coordinates and the basic topological information.
    
    The quantized coordinates of the vertex with index {uxv} will be {vpos[uxv]}. 
    
    The endpoints of the unoriented edge with index {uxe} will be
    the vertices with indices {endv[uxe].c[0..1]}.  
    
    The sides of the unoriented face with index {uxf} will be the 
    unoriented edges with indices {side[uxf].c[0..2]}.  
    
    The procedure also computes derived informetion
    such as edge degrees, the fields {.minZ,.maxZ} of each triangle,
    and the bounding box of the mesh.
    
    If {checkSorted} is true, the procedure verifies 
    that the faces are sorted in order of non-decreasing
    {.minZ} field. */

/* SIZE LIMITS */

#define stmesh_n_MAX ((uint32_t)(1u << 30))
  /* Maximum number of vertices, oriented edges, or oriented faces in a
    mesh. An index of any of those things will safely fit in an {int32_t}
    as well as in an {uint32_t}. */

#define stmesh_nv_MAX (stmesh_n_MAX)
  /* Max number of vertices in a mesh. */

#define stmesh_ne_MAX (stmesh_n_MAX / 2)
  /* Max number of UNORIENTED edges in a mesh. */

#define stmesh_nf_MAX (stmesh_n_MAX / 6)
  /* Max number of UNORIENTED triangles in a mesh. */

#endif
