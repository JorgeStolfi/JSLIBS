/* Arbitrary polygonal figure, with quasimanifold topology. */
/* Last edited on 2016-04-21 18:40:07 by stolfilocal */

#ifndef stpoly_H
#define stpoly_H

#define _GNU_SOURCE
#include <stdio.h>

#include <vec.h>
#include <bool.h>
#include <i2.h>
#include <r2.h>

#include <stpoly_STP.h>

typedef struct stpoly_rep_t *stpoly_t;
  /* An {stpoly_t} is a pointer to a data structure that represents a
    read-only list of segments (/edegs/) in {\RR^2}, with topological information.
    It has a list of distinct vertices and a list of distinct edges. 
    The coordinates of each vertex are rounded to even
    multiples of a fundamental length unit {eps}, stored in the
    structure. */
     
float stpoly_get_eps(stpoly_t poly);
  /* Returns the fundamental length unit (in mm) of the given polygonal figure. */
    
void stpoly_get_bounding_box(stpoly_t poly, i2_t *minQP, i2_t *maxQP);
  /* Stores into {*minQP,*maxQP} the lower and upper
    corner of the bounding box of {poly}. */

void stpoly_print_bounding_box(FILE *wr, stpoly_t poly);
  /* Writes to {wr} a readable display of the poly's bounding box,
    in quantized coordinates and in {mm}. */

uint32_t stpoly_vert_count(stpoly_t poly);
uint32_t stpoly_edge_count(stpoly_t poly);
  /* These procedures return the number of distinct vertices
    and edges in the given figure {poly}, respectively. */

/* VERTICES AND EDGES. */

/* The parts of a polygonal figure - vertices and edges -- are represented
  by records with private format. The edge pointer types below
  may have their low-order bits modified to indicate an orientation for
  the element. Thus they must not be dereferenced directly. */

typedef struct stpoly_vert_rep_t *stpoly_vert_t;
  /* A pointer to a record that represents a vertex in a {stpoly_t}. */

typedef struct stpoly_edge_rep_t *stpoly_edge_t;
  /* A pointer to a record that represents an ORIENTED edge in a {stpoly_t};
    that is, an edge with a specific sense of traversal.  The orientation defines
    which endpoint is the /origin/ and which is the /destination/. */

/* VERTEX OPERATIONS */    
 
i2_t stpoly_vert_get_pos(stpoly_vert_t v);
  /* Returns the quantized coordinates of the vertex {v}. */

r2_t stpoly_unround_point(i2_t *p, float eps);
  /* Converts the quantized point {*p} to a non-quantized point,
    with coordinates in mm, by multiplying it by {eps}. */
    
uint32_t stpoly_vert_degree(stpoly_vert_t v);
  /* Returns the degree of vertex {v}, that is, the number of 
    edges incident to {v}. */

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

void stpoly_edge_get_yrange(stpoly_edge_t e, int32_t *minYP, int32_t *maxYP);
  /* Returns in {*minYP,*maxYP} the minimum and maximum quantized {Y} coordinates of the 
    endpoits of the edge {e}.  Ignores the edge's orientation. */
 
/* GEOMETRIC SLICING OPERATIONS */

bool_t stpoly_edge_crosses_line(stpoly_edge_t e, int32_t pP);
  /* Returns {TRUE} iff the edge {e} crosses the horizontal line 
    with quantized {Y}-coordinate {pY}. The outcome is undefined 
    if the line contains one or both endpoints of the edge. */
 
double stpoly_edge_line_intersection(stpoly_edge_t e, int32_t pY, double eps);
  /* Computes the coordinates (in mm) of the point where the edge {e}
    intersects the horizontal line with quantized {Y}-coordinate {pY}.
    The outcome is undefined if the line does not intersect the edge,
    or contains one or both endpoints of the edge. */

/* ELEMENT INDEXING */

typedef uint32_t stpoly_vert_unx_t;
  /* An integer that identifies a
    vertex of a poly, namely the index of its vertex record in 
    in the table {poly.v}.  Must be in {0..nv-1}, where 
    {nv=stpoly_vert_count(poly)}. */

typedef uint32_t stpoly_edge_unx_t;
  /* An integer that identifies an UNORIENTED edge in a
    {stpoly_t}, namely the index of its vertex record in in the table
    {poly.e}. Must be in the range {0..ne-1}, where
    {ne=stpoly_edge_count(poly)}. */

stpoly_vert_t stpoly_get_vert(stpoly_t poly, stpoly_vert_unx_t uxv);
stpoly_edge_t stpoly_get_edge(stpoly_t poly, stpoly_edge_unx_t uxe);
  /* Get the vertex or edge of {poly} with the speficied index.
    The index must be in the appropriate range. Edges are
    returned in their natural orientation. */

stpoly_vert_unx_t stpoly_vert_get_unx(stpoly_t poly, stpoly_vert_t v);
stpoly_edge_unx_t stpoly_edge_get_unx(stpoly_t poly, stpoly_edge_t e);
  /* Returns the index of the vertex, edge, or face in {poly}.
    disregards the orientation of edges. */

/* CONSTRUCTION */

stpoly_t stpoly_new_desc(float eps, uint32_t nv_max, uint32_t ne_max);
  /* Creates a new descriptor record for a poly, with capacity for
    {nv_max} vertices and {ne_max} edges,
    but with no vertices and no edges. */

stpoly_vert_unx_t stpoly_add_vert(stpoly_t poly, i2_t *pos);
  /* Adds a vertex to {poly}, given its quantized coordinates {*pos}.
    Returns its index in the poly. Updates the bounding box of the poly
    to include that vertex. The vertex had better be used later in some
    edge. */

stpoly_edge_unx_t stpoly_add_edge(stpoly_t poly, stpoly_vert_unx_t uxv[]);
  /* Adds an edge to the poly, given the indices {uxv[0..1]} of its
    endpoint vertices. The indices must be in increasing order. Sets its
    degree to zero. Returns the index of the edge in the poly. The edge
    had better be used later by some face. Assumes that {poly->v} is set. */

void stpoly_free(stpoly_t poly);
  /* Releases all heap storage allocated for {poly}, including the 
    descriptor record. */

/* MISCELLANEOUS */

void stpoly_print_vert_degrees(FILE *wr, stpoly_t poly);
  /* Counts how many vertices have degree {0,1,...}, 
    prints them to {wr}.  Aborts with error message 
    if there is a vertex with degree 0. */ 

void stpoly_print(FILE *wr, stpoly_t poly);
void stpoly_vert_print(FILE *wr, stpoly_t poly, stpoly_vert_t v);
void stpoly_edge_print(FILE *wr, stpoly_t poly, stpoly_edge_t e);
  /* Writes legibly to {wr} the address, number, and fields of the element. */

void stpoly_check(stpoly_t poly);
void stpoly_vert_check(stpoly_t poly, stpoly_vert_t v);
void stpoly_edge_check(stpoly_t poly, stpoly_edge_t e);
  /* Runs some consistency checks on the topological operations
    on those elements. */

/* INPUT/OUTPUT */

stpoly_t stpoly_read_STP(char *fileName, bool_t binary, float eps, uint32_t neGuess, bool_t even);
  /* Reads the STP file with given name (ascii or binary)
    and converts it to a {stpoly_t} structure. 

    See {stpoly_read_STP_INFO} below for an explanation of how
    the topology of the poly is recovered from the unstructured STP file.

    The {neGuess} parameter is a hint for the number of (unoriented)
    edges in the figure. It is used to pre-allocate tables of edges and
    other items. It can be any number, even zero; however, the procedure
    is more efficient if {nf_guess} is equal to the actual number of edges, or
    slightly higher.
    
    If {even} is true, each vertex coordinate is quantized by rounding
    to the nearest *even* multiple of the fundamental length {eps}. If {even} is
    false, it is rounded to the nearest integer multiple of {eps}, even or odd. */
    
#define stpoly_read_STP_INFO \
  "The topology is determined by first rounding every coordinate of every vertex (segment endpoint) to" \
  " an integer multiple of some length unit {EPS} specified by the user.  Any" \
  " endpoints that become coincident by such rounding are assumed to be the same vertex" \
  " of the polygonal figure.  Any edge (segment) that has two coincident endpoints is flagged" \
  " and discarded.\n" \
  "\n" \
  "  Every vertex" \
  " will be incident to at least one edge, usually two, possibly" \
  " three or more. It is a fatal error if there are two edges" \
  " with exactly the same vertices after quantization.  Reducing the length unit {EPS} may get" \
  " around such problem; otherwise, the input STP file should be edited," \
  " by removing *both* ofending edges."

/* SIZE LIMITS */

#define stpoly_n_MAX ((uint32_t)(1u << 30))
  /* Maximum number of vertices  or oriented edges in a
    polygonal figure. An index of any of those things will safely fit in an {int32_t}
    as well as in an {uint32_t}. */

#define stpoly_nv_MAX (stpoly_n_MAX)
  /* Max number of vertices in a poly. */

#define stpoly_ne_MAX (stpoly_n_MAX / 2)
  /* Max number of UNORIENTED edges in a poly. */

#endif
