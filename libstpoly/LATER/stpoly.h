/* Polygonal regions. */
/* Last edited on 2016-04-13 11:27:04 by stolfilocal */

#ifndef stpoly_H
#define stpoly_H

#define _GNU_SOURCE
#include <stdio.h>

#include <vec.h>
#include <bool.h>
#include <i3.h>
#include <r2.h>
#include <r3.h>

#include <stpoly_STP.h>

/* TRIANGLE MESHES WITH TOPOLOGICAL INFORMATION */

typedef struct stpoly_rep_t *stpoly_t;
  /* An {stpoly_t} is a pointer to a data structure that represents a
    read-only list of triangles in {\RR^3}, with topological
    information. It has a list of distinct vertices, a list of distinct
    edges, and a list of distinct faces (triangles). The coordinates of
    each vertex are rounded to even multiples of a fundamental length
    unit {eps}, stored in the structure. */
     
float stpoly_get_eps(stpoly_t mesh);
  /* Returns the fundamental length unit (in mm) of the given mesh. */
    
void stpoly_get_bounding_box(stpoly_t mesh, i3_t *minQP, i3_t *maxQP);
  /* Stores into {*minQP,*maxQP} the lower and upper
    corner of the bounding box of {mesh}. */

void stpoly_print_bounding_box(FILE *wr, stpoly_t mesh);
  /* Writes to {wr} a readable display of the mesh's bounding box,
    in quantized coordinates and in {mm}. */

uint32_t stpoly_vert_count(stpoly_t mesh);
uint32_t stpoly_edge_count(stpoly_t mesh);
uint32_t stpoly_face_count(stpoly_t mesh);
  /* These procedures return the number of distinct vertices,
    edges, and faces in the given {mesh}, respectively. */

/* VERTICES, EDGES, AND FACES. */

/* The parts of the mesh - vertices, edges, and faces -- are represented
  by records with private format. The edge and face pointer types below
  may have their low-order bits modified to indicate an orientation for
  the element. Thus the may not be dereferenced directly. */

typedef struct stpoly_vert_rep_t *stpoly_vert_t;
  /* A pointer to a record that represents a vertex in a {stpoly_t}. */

typedef struct stpoly_edge_rep_t *stpoly_edge_t;
  /* A pointer to a record that represents an ORIENTED edge in a {stpoly_t};
    that is, an edge with a specific sense of traversal.  The orientation defines
    which endpoint is the /origin/ and which is the /destination/,
    and also what is the circular ordering of the faces around the edge, if there are
    more than two (by the right-hand rule). */
     
typedef struct stpoly_face_rep_t *stpoly_face_t;
  /* A pointer to a record that represents an ORIENTED triangle in a
    {stpoly_t}. The orientation specifies one of the two possible
    senses of traversal of its perimeter, and also designates a specific
    side as being the triangle's /base edge/. Thus, there are six distincg
    oriented triangles for each actual (unoriented) triangle. */

/* VERTEX OPERATIONS */    
 
i3_t stpoly_vert_get_pos(stpoly_vert_t v);
  /* Returns the quantized coordinates of the vertex {v}. */

r3_t stpoly_unround_point(i3_t *p, float eps);
  /* Converts the quantized point {*p} to a non-quantized point,
    with coordinates in mm, by multiplying it by {eps}. */

/* EDGE OPERATIONS */

stpoly_edge_t stpoly_edge_reverse(stpoly_edge_t e, int k);
  /* If {k} is odd, returns the oriented edge which is the same unoriented edge as
    {e}, but with the opposite orientation.  If {k} is even, retruns {e}
    itself. */

stpoly_edge_t stpoly_edge_natural(stpoly_edge_t e);
  /* Returns the same edge as {e}, in its "natural" orientation.
    That is, {stpoly_edge_natural(e)==stpoly_edge_natural(e')=e*}
    for any oriented edge {e}, where {e'} is the reverse of {e},
    and {e*} is either {e} or {e'}. */

stpoly_vert_t stpoly_edge_get_endpoint(stpoly_edge_t e, int k);
  /* Returns the vertex that is the origin ({k = 0}) or
    destination ({k = 1}) of {e}, taking its orientation into account. */

void stpoly_edge_get_endpoints(stpoly_edge_t e, stpoly_vert_t v[]);
  /* Stores into {v[0..1]} the end vertices of {e}, in the order determined
    by its orientation.  Namely {v[0]}  will be the origin, and {v[1]} the destintion. */
    
uint32_t stpoly_edge_degree(stpoly_edge_t e);
  /* Returns the degree of edge {e}, that is, the number of 
    faces incident to {e}. Ignores the orientation. */

/* FACE OPERATIONS */    

stpoly_edge_t stpoly_face_get_base(stpoly_face_t f);
  /* Return the oriented base edge of the oriented face {f}. */
  
stpoly_face_t stpoly_face_flip(stpoly_face_t f, int k);
  /* If {k} is odd, returns the oriented face which is the same
    unoriented face as {f}, but taken with the opposite orientation. If
    {k} is even, returns {f} itself.
    
    When the orientation is flipped, the base edge of the result will be
    the base edge of {f}, reversed. The new orientation specifies that
    the sides of {f} will be traversed in the opposite order, starting from
    that new base edge, and each side will traversed in the opposite
    sense. */

stpoly_face_t stpoly_face_shift(stpoly_face_t f, int k);
  /* Returns the oriented face that is the same unoriented
    face as {f}, with the perimeter traversed in the same sense, but whose
    base edge is {k} sides ahead of {f}'s base edge, in that sense of traversal.
    
    In particular, {stpoly_face_shift(f,k)} is the same as
    {stpoly_face_shift(f,k+3)}, and
    {stpoly_face_shift(stpoly_face_flip(f),k)} is
    {stpoly_face_flip(stpoly_face_shift(f,-k))}, for all {k}. */

stpoly_face_t stpoly_face_natural(stpoly_face_t f);
  /* Returns the oriented face that is the same unoriented
    face as {f}, in its "natural" sense and base edge.
    That is, {stpoly_face_natural(f)==stpoly_face_natural(f')=f*}
    for any oriented versions {f} and {f'} of the same 
    unoriented face; where {f*} is a particular
    flipped and/or shifted version of {f}. */

void stpoly_face_get_sides(stpoly_face_t f, stpoly_edge_t e[]);
  /* Returns oriented edges that are the sides of the 
    face {f}, starting with its base edge and progressing along the perimeter in the
    direction specified by its orientation.  The edges will be oriented so as to 
    agree with that sense of traversal. */

void stpoly_face_get_corners(stpoly_face_t f, stpoly_vert_t v[]);
  /* Stores in {v[0..2]} the three corners of the oriented face {f},
    starting with the vertex opposte to its base edge and progressing along the perimeter in the
    direction specified by its orientation. */

void stpoly_face_get_zrange(stpoly_face_t f, int32_t *minZP, int32_t *maxZP);
  /* Returns in {*minZP,maxZP} the minimum and maximum quantized {Z} coordinates of the 
    corners of the triangle.  Ignores the triangle's orientation. */
 
/* GEOMETRIC SLICING OPERATIONS */

bool_t stpoly_edge_crosses_plane(stpoly_edge_t e, int32_t pZ);
  /* Returns {TRUE} iff the edge {e} crosses the horizontal plane 
    with quantized {Z}-coordinate {pZ}. The outcome is undefined 
    if the plane contains one or both endpoints of the edge. */
 
r2_t stpoly_edge_plane_intersection(stpoly_edge_t e, int32_t pZ, double eps);
  /* Computes the coordinates (in mm) of the point where the edge {e}
    intersects the slicing plane with quantized {Z}-coordinate {pZ}.
    The outcome is undefined if the plane does not intersect the edge,
    or contains one or both endpoints of the edge. */

void stpoly_face_get_sliced_sides(stpoly_t mesh, stpoly_face_t f, int32_t pZ, stpoly_edge_t e[]);
  /* Given the quantized {Z}-coordinate {pZ} of a slicing plane, finds
    the two sides of the oriented face {f} that are crossed by the
    plane. Assumes that the plane intersects the face but does not go
    through any of its vertices. Returns the two edges in arbitrary order,
    in their natural orientation. */

/* ELEMENT INDEXING */

typedef uint32_t stpoly_vert_unx_t;
  /* An integer that identifies a
    vertex of a mesh, namely the index of its vertex record in 
    in the table {mesh.v}.  Must be in {0..nv-1}, where 
    {nv=stpoly_vert_count(mesh)}. */

typedef uint32_t stpoly_edge_unx_t;
  /* An integer that identifies an UNORIENTED edge in a
    {stpoly_t}, namely the index of its vertex record in in the table
    {mesh.e}. Must be in the range {0..ne-1}, where
    {ne=stpoly_edge_count(mesh)}. */

typedef uint32_t stpoly_face_unx_t;
  /* An integer that identifies an UNORIENTED triangle in a
    {stpoly_t}, namely the index of its vertex record in in the table
    {mesh.f}. Must be in the range {0..nf-1}, where {ne =
    stpoly_edge_count(mesh)}. */

stpoly_vert_t stpoly_get_vert(stpoly_t mesh, stpoly_vert_unx_t uxv);
stpoly_edge_t stpoly_get_edge(stpoly_t mesh, stpoly_edge_unx_t uxe);
stpoly_face_t stpoly_get_face(stpoly_t mesh, stpoly_face_unx_t uxf);
  /* Get the vertex, edge, or face of {mesh} with the speficied index.
    The index must be in the appropriate range. Edges and faces are
    returned in their natural orientation. */

stpoly_vert_unx_t stpoly_vert_get_unx(stpoly_t mesh, stpoly_vert_t v);
stpoly_edge_unx_t stpoly_edge_get_unx(stpoly_t mesh, stpoly_edge_t e);
stpoly_face_unx_t stpoly_face_get_unx(stpoly_t mesh, stpoly_face_t f);
  /* Returns the index of the vertex, edge, or face in {mesh}.
    disregards the orientation of edges and faces. */

/* CONSTRUCTION */

stpoly_t stpoly_new_desc(float eps, uint32_t nv_max, uint32_t ne_max, uint32_t nf_max);
  /* Creates a new descriptor record for a mesh, with capacity for
    {nv_max} vertices, {ne_max} edges, and {nf_max} faces,
    but with no edges, no vertices, and no faces. */

stpoly_vert_unx_t stpoly_add_vert(stpoly_t mesh, i3_t *pos);
  /* Adds a vertex to {mesh}, given its quantized coordinates {*pos}.
    Returns its index in the mesh. Updates the bounding box of the mesh
    to include that vertex. The vertex had better be used later in some
    edge. */

stpoly_edge_unx_t stpoly_add_edge(stpoly_t mesh, stpoly_vert_unx_t uxv[]);
  /* Adds an edge to the mesh, given the indices {uxv[0..1]} of its
    endpoint vertices. The indices must be in increasing order. Sets its
    degree to zero. Returns the index of the edge in the mesh. The edge
    had better be used later by some face. Assumes that {mesh->v} is set. */
        
stpoly_face_unx_t stpoly_add_face(stpoly_t mesh, stpoly_edge_unx_t uxe[]);
  /* Adds a face to the mesh, given the indices {uxe[0..2]} of the
    unoriented edges that bound the face. The indices should be sorted
    in increasing order. Increments the degree of the edges that are
    sides of {f}. Returns the index of the face in the mesh. Assumes
    that {mesh.{nv,v,ne,e}} are set. */

/* MISCELLANEOUS */

void stpoly_print_edge_degrees(FILE *wr, stpoly_t mesh);
  /* Counts how many edges have degree {0,1,...}, 
    prints them to {wr}.  Aborts with error message 
    if there is an edge with degree 0. */ 

void stpoly_print(FILE *wr, stpoly_t mesh);
void stpoly_vert_print(FILE *wr, stpoly_t mesh, stpoly_vert_t v);
void stpoly_edge_print(FILE *wr, stpoly_t mesh, stpoly_edge_t e);
void stpoly_face_print(FILE *wr, stpoly_t mesh, stpoly_face_t f);
  /* Writes legibly to {wr} the address, number, and fields of the element. */

void stpoly_check(stpoly_t mesh);
void stpoly_vert_check(stpoly_t mesh, stpoly_vert_t v);
void stpoly_edge_check(stpoly_t mesh, stpoly_edge_t e);
void stpoly_face_check(stpoly_t mesh, stpoly_face_t f);
  /* Runs some consistency checks on the topological operations
    on those elements. */

/* INPUT/OUTPUT */

stpoly_t stpoly_read_STP(char *fileName, bool_t binary, float eps, uint32_t nf_guess, bool_t even, bool_t checkSorted);
  /* Reads the STP file with given name (ascii or binary)
    and converts it to a {stpoly_t} structure. 

    See {stpoly_read_STP_INFO} below for an explanation of how
    the topology of the mesh is recovered from the unstructured STP file.

    The {nf_guess} parameter is a hint for the number of (unoriented)
    faces in the mesh. It is used to pre-allocate tables of faces and
    other items. It can be any number, even zero; however, the procedure
    is more efficient if {nf_guess} is equal to the number of faces, or
    slightly higher.
    
    If {even} is true, each vertex coordinate is quantized by rounding
    to the nearest *even* multiple of the fundamental length {eps}. If {even} is
    false, it is rounded to the nearest integer multiple of {eps}, even or odd.

    if {checkSorted} is true, checks whether the triangles in {mesh.f}
    are sorted by their {.minZ} field in non-decreasing order. Aborts
    with message if not. */
    
#define stpoly_read_STP_INFO \
  "The topology is determined by first rounding every coordinate of every vertex to" \
  " an integer multiple of some langth unit {EPS} specified by the user.  Triangle" \
  " corners that become coincident by such rounding are assumed to be the same vertex" \
  " of the mesh.  Any triangle that has one or two coincident vertices is flagged" \
  " and discarded.  It is a fatal error if two triangles with exactly the same vertices.\n" \
  "\n" \
  "  Then any two sides of two different triangles that have the same" \
  " vertices as endpoints are assumed to be the same edge of the mesh.  Every edge" \
  " then must be incident to at least one face, usually two, possibly three or more.  Every" \
  " vertex too must be incident to at least two distinct edges and at least one triangle.\n" \
  "\n" \
  "  It is a fatal error if two triangles of the mesh have the same three" \
  " vertices after quantization. Reducing the length unit {EPS} may get" \
  " around such problem; otherwise, the input STP file should be edited" \
  " by removing *both* ofending triangles."

/* SIZE LIMITS */

#define stpoly_n_MAX ((uint32_t)(1u << 30))
  /* Maximum number of vertices, oriented edges, or oriented faces in a
    mesh. An index of any of those things will safely fit in an {int32_t}
    as well as in an {uint32_t}. */

#define stpoly_nv_MAX (stpoly_n_MAX)
  /* Max number of vertices in a mesh. */

#define stpoly_ne_MAX (stpoly_n_MAX / 2)
  /* Max number of UNORIENTED edges in a mesh. */

#define stpoly_nf_MAX (stpoly_n_MAX / 6)
  /* Max number of UNORIENTED triangles in a mesh. */

#endif
