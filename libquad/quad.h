#ifndef quad_H
#define quad_H

/* The quad-edge data structure (oriented surface version). */
/* Last edited on 2024-12-05 10:39:45 by stolfi */

#define quad_H_copyright \
  "Copyright © 1996, 2006 State University of Campinas (UNICAMP).\n\n" jslibs_copyright "\n\n" \
  "NOTE: this copyright notice does not claim to supersede any copyrights" \
  " that may apply to the original DEC implementation of the quad-edge" \
  " data structure."

/* The quad-edge data structure encodes the topology of a graph drawn on
  an orientable compact 2-D manifold (i.e. a borderless surface of
  finite extent) in such a way that every face is a topological disk.
  For details, see
 
    "Primitives for the Manipulation of General Subdivisions 
    and the Computation of Voronoi Diagrams"
 
    L. Guibas, J. Stolfi, ACM Transcations on Graphics, April 1985
 
  Originally written by Jim Roth (DEC CADM Advanced Group) in May/1986.
  Modified by J. Stolfi in Apr/1993. Changed by J. Stolfi in Dec/2011 to
  be procedure-based instead of macro-based. */

#include <stdio.h>
#include <stdint.h>

#include <jslibs_copyright.h>
#include <bool.h>
#include <vec.h>

typedef struct quad_rep_t *quad_arc_t; 
  /* A directed edge, primal or dual.  A null arc is {NULL} */

quad_arc_t quad_nop(quad_arc_t e);
  /* The identity or `no operation' function that returns {e}
    itself; it is provided here for completeness. */ 

quad_arc_t quad_sym(quad_arc_t e);
  /* The {sym} operator reverses the longitudinal and transversal
    orientations. Intuitively, it rotates {e} by a half-turn. Thus
    {sym(e) != e} but {sym(sym(e)) = e}, for any {e}. */
     
/* DUALITY OPERATORS

  Actually, at any moment the quad-edge data structure represents
  *two* maps on {X}, each the dual of the other.
  
  The two maps are disjoints, and there is a one-to-one /duality/
  correspondence {*} between their elements such that, for any element
  {A} of either map, the dimensions of {A} and {A*} add to 2.
  Moreover, {A*} intersects {A} at precisely one point, and is
  disjoint from {B} or {B*}, for any other element {B!= A} with the
  same dimension as {A}. */

quad_arc_t quad_rot(quad_arc_t e);
quad_arc_t quad_tor(quad_arc_t e);
  /* The {rot} operation returns an arc {f} on the dual edge {F} of
    the edge {ed} underlying {e}, oriented so that {f} crosses {e}
    leftwards and {e} crosses {f} rightwards. Intuitively, {rot}
    rotates {e} by a positive quarter-turn. Thus, {rot(e) != e},
    {rot(e) != sym(e)}, but {rot(rot(e)) == sym(e)}, for any {e}.
    
    The {tor} operator is the inverse of {rot}: it rotates {e} by a
    negative quarter-turn. It is equivalent to {rot(sym(e))}.
    
    Both {rot} and {tor} preserve the circular orientation 
    of the arc itself. */
  
/* TUMBLING OPERATORS 
 
  The /tumbling operators/ are the four functions {nop}, {rot},
  {sym}, {tor}. This set is closed under composition: any 
  sequence of tumblings is equivalent to  a single tumbling. 
  Namely, the composition {G(F(e)) = e.F.G} is given by the table

            | G                                           
            | nop   sym    rot   tor   
    --------|--------------------------
    F nop   | nop   sym    rot   tor   
      sym   | sym   nop    tor   rot   
            |                    
      rot   | rot   tor    sym   nop    
      tor   | tor   rot    nop   sym   

  Note that the operators {sym,rot,tor} commute. */

quad_arc_t quad_rot_n(quad_arc_t e, int32_t r);
  /* This generic operator applies {r} times the {quad_rot} to arc {e}. */

/* VERTEX/FACE WALKING OPERATORS

  Each of the operators in this section steps through a circular list
  {S} of arcs that includes the given arc {e}. All arcs in the list
  {S} have a common /pivot element/: the same origin vertex, the same
  destination vertex, the same left face, or the same right face.

  The order of enumeration of {S} is defined by a circular
  orientation {C} on the pivot element (that is, on the pivot face,
  or on a neighborhood of the pivot vertex), which is derived from
  the circular orientation of {e}.

  The first letter of the operator's name specifies the pivot element
  ('{o}' for origin, '{d}' for destination, '{l}' for left face, and '{r}'
  for right face). It also specifies how the circular orientation of
  {e} is transferred to the pivot: by extrapolation through the
  origin ('{o}') or destination ('{d}') extremity, or into the left ('{l}')
  or right ('{r}') bank of the arc.

  Note that the same edge may appear twice in the list {S}. In such
  cases, the two occurences will have opposite orientations. On the other hand, in
  ordinary quad-edge structures the list {S} of {e} will never include
  any of its dual versions ({rot(e)} or {tor(e)}). */

quad_arc_t quad_onext(quad_arc_t e);
quad_arc_t quad_oprev(quad_arc_t e);
 /* The next and previous arc, respectively,
   among the arcs with same origin as {e}.  */

quad_arc_t quad_dnext(quad_arc_t e);
quad_arc_t quad_dprev(quad_arc_t e);
 /* The next and previous arc, respectively,
   among the arcs with same destination as {e}.  */

quad_arc_t quad_lnext(quad_arc_t e);
quad_arc_t quad_lprev(quad_arc_t e);
 /* The next and previous arc, respectively,
   among the arcs with same left face as {e}.  */

quad_arc_t quad_rnext(quad_arc_t e);
quad_arc_t quad_rprev(quad_arc_t e);
 /* The next and previous arc, respectively,
   among the arcs with same right face as {e}.  */

quad_arc_t quad_walk(quad_arc_t e, int32_t r, int32_t n);
  /* Returns {tor^r(onext^n(rot^r(e)))}.
    For negative values of {r} and/or {f}, the result is defined by
    the identities {walk(e,-r,n) = rot^r(onext^n(tor^r(e)))} and
    {walk(e,r,-n) = tor^r(oprev^n(rot^r(e)))}. Thus, for example, 
    
      {walk(e,0,1) == onext(e)}, 
      {walk(e,1,1) == rnext(e)}, 
      {walk(e,2,1) == dnext(e)}, 
      {walk(e,3,1) == lnext(e)}, 
      {walk(e,0,-1) == oprev(e)}, 
      {walk(e,2,3) == dnext(dnext(dnext(e)))}, 
      {walk(e,-5,-2) == lprev(lprev(e))}, 
      
    et coetera. */
    
/* Topology_changing operators: */

quad_arc_t quad_make_edge(void);
  /* Creates a structure consisting of an isolated non-loop primal edge
    and its dual (which is a loop). Returns one of the orientations
    of the primal edge. */

void quad_destroy_edge(quad_arc_t e);
  /* Releases the edge record of the arc {e} and of its 
    dual, in all oriented variants. */

void quad_splice(quad_arc_t a, quad_arc_t b);
  /* Applies the {splice} operator to the directed 
    arcs {a} and {b}.  The result is undefined if one
    arc is reachable from the other by a sequence 
    of {rot} and {onext} operations with an odd
    number of {rot}s. */

/* EDGE RECORDS */

typedef struct quad_edge_rec_t *quad_edge_t;
  /* An edge record, representing an undirected and unoriented edge {ed}
    of the primal map together with the corresponding edge {F} of
    the dual map. It is shared by the four arcs that consist of {ed} or {F}
    taken with both logitudinal orientations. */

quad_edge_t quad_edge(quad_arc_t e);
  /* Obtains the edge record of an arc reference {e}. Satisfies
    {edge(rot(e)) == edge(e)}.  */
  
void *quad_odata(quad_arc_t e);
void *quad_rdata(quad_arc_t e);
void *quad_ddata(quad_arc_t e);
void *quad_ldata(quad_arc_t e);
  /* Each directed arc has an associated "data" field (an arbitrary pointer)
    which can be thought of as being stored near its origin end.
    These procedures return the data field associated with {e},
    {rot(e)}, {sym(e)}, and {tor(e)}, respectively.  Note that 
    {odata(e)} and {ddata(e)} are distinct variables even when the
    arc is a loop. */

void quad_set_odata(quad_arc_t e, void *p);
void quad_set_rdata(quad_arc_t e, void *p);
void quad_set_ddata(quad_arc_t e, void *p);
void quad_set_ldata(quad_arc_t e, void *p);
  /* These procedures store {p} into the data value associated with {e},
    {rot(e)}, {sym(e)}, and {tor(e)}, respectively. */

/* TUMBLE CODE

  A value {e} of type {arc_t} consists of a pointer
  {ed = edge(e)} to an /edge record/, of type {edge_rec_t}, 
  and a two-bit /tumble code/ {tc = tumble_code(e)}.
  
  The pointer {ed} and the tumble code {tc} are packed together as a
  single {void *} value. This trick assumes that the address of any
  edge record is a multiple of 4 bytes (32 bits), so {tc} can be stored
  in its lowest two bits, without ambiguity. Thus, although an {arc_t}
  is formally an address, it should not be used as such, since it
  generally points to some random byte inside the edge record. */

typedef uint8_t quad_bits_t;
   /* A data type used to hold a few bits (usually at the low-order end). */
    
quad_bits_t quad_tumble_code(quad_arc_t e);
  /* Returns a two-bit code (an integer in {0..3}) that identifies
    the arc {e} among the four arcs with the same edge record. Two
    arcs {e,f} are identical if and only if {edge(e)==edge(f)} and
    {tumble_code(e)==tumble_code(f)}. */
  
quad_arc_t quad_orient(quad_edge_t ed, quad_bits_t tc);
  /* Returns the unique arc that has edge record {ed} and tumble code {tc}. 
    Only the two lower-order bits of {tc} are used.
    
    This procedure can be convenient in special situations, like
    reading or writing a quad-edge structure, or enumerating all arcs
    with the same edge. It is not very useful in general because the
    meaning of the tumble codes is implementation-dependent. In
    particular, {orient(edge(e),tumble_code(f))} is hardly sensible
    if {edge(f)!=edge(e)}. */

/* ORIENTATION BITS
  
  The tumble code of an arc {e} can be mapped to two Boolean
  attributes: the /longitudinal orientation bit/ {lon_bit(e)}, the
  /dual-primal bit/ {dop_bit(e)}. These three bits completely determine the tumble
  code (in some implementation-dependent way). */
  
quad_bits_t quad_lon_bit(quad_arc_t e);
  /*  The /longitudinal orientation bit/ of {e}. This attribute is
    reversed by {sym}, and may or may not be changed by {rot}. */
  
quad_bits_t quad_dop_bit(quad_arc_t e); 
  /* The /dual-primal bit/ of {e}. This attribute is
    reversed by {rot} and preserved by {sym}. */

/* CAVEAT ON ORIENTATION BITS 

  It is important to note that the tumbing code and the above bits do
  not distinguish the orientations or the primal/dual character in a
  global sense, but only among the arcs of the same pair of mutually
  dual edges. Thus, for example, the arcs {e} and {onext(e)} are both
  in the primal map, or both in the dual map; but their underlying
  edges are (generally) distinct, so they may have different
  {dop_bit}s. In general, there is no constraint on the orientation
  bits of two arcs on distinct edges, even when they are connected by
  a walking function such as {onext}, {rnext}, etc.
  
  There is also no constraint on the {sym_bit}s of {e} and {rot(e)}.
  In fact, if {sym_bit(e)==sym_bit(rot(e))} for some arc {e}, then one
  must have {sym_bit(f)!=sym_bit(rot(f))} for the arc {f==rot(e)}. */

/* NULL ARCS */

#define quad_arc_NULL ((quad_arc_t)NULL)
  /* A special {arc_t} value that means `no such arc'. */

bool_t quad_arc_is_null(quad_arc_t e);
  /* TRUE if {e} is a null {quad_arc_t}, equivalent to
    {quad_edge(e)==NULL}. Note that this is a stronger test than {e ==
    NULL} or {e == quad_arc_NULL} because one may have tumbled
    versions of the {NULL} edge. */

/* VECTORS OF ARCS AND EDGES */

vec_typedef(quad_arc_vec_t,quad_arc_vec,quad_arc_t);
  /* An {quad_arc_vec_t} is a vector of {quad_arc_t}s. */

vec_typedef(quad_edge_vec_t,quad_edge_vec,quad_edge_t);
  /* A {quad_edge_vec_t} is a vector of {quad_edge_t}s. */

/* Map traversal: */

void quad_enum(quad_arc_vec_t *root, void visit_proc(quad_arc_t e));
  /* Enumerates undirected primal edges reachable from {root.e[0..root.ne-1]}.
    
    Calls visit_proc(e, closure) for every arc {e} that can be reached from
    a {root} arc by a chain of {SYM} and {ONEXT} calls; except that exactly one
    of {e} and {SYM(e)} is visited. */

/* Edge numbers: */

typedef uint64_t quad_edge_id_t; 
  /* An integer stored in an edge record. */

quad_edge_id_t quad_edge_id(quad_edge_t ed);
  /* Returns the current number of edge {ed}. */

void quad_set_edge_id(quad_edge_t ed, quad_edge_id_t eid);
  /* Sets the number of edge {ed} to {eid}. */

quad_edge_id_t quad_renumber_edges(quad_arc_vec_t *root, quad_arc_vec_t *A);
  /* Renumbers all edge records reachable from {root[0..nr-1]}
    sequentially from 0, where {nr = root->ne}. Returns the number {ne}
    of edges found. 
    
    If {A} is not NULL, stores into each element {A->e[k]} with {k} in
    {0..ne} one reachable arc from the edge {ed} such that
    {quad_edge_id(ed) == k}. The vector {A} is expanded as needed and
    trimmed to {ne} elements. */

/* ARC INPUT/OUTPUT */

void quad_write_arc(FILE *wr, quad_arc_t e, int32_t width);
  /* Writes the arc {e} to {wr} in the format "{eid}:{tc}" where {eid}
    is {quad_edge_id(edge(e))}, and {tc} is {quad_tumble_code(e)} in binary.
    The {eid} is padded with spaces to have at least {width}
    characters. */

quad_arc_t quad_read_arc(FILE *rd, quad_arc_vec_t *A);
  /* Parses from file {rd} a description of an arc, in the format
    "{eid}:{tc}" where {eid} is an edge number and {tc} is a tumble
    code. Requires an edge table {A}, which is a vector of arcs such
    that {A->e[i]} is any arc on an edge edge {ed} with {quad_edge_id(ed) == i}.
    The resulting arc {e} has {quad_edge(e) == quad_edge(A->e[eid])} and
    {quad_tumble_code(e) == tc}. */

/* MAP INPUT/OUTPUT */

void quad_write_map(FILE *wr, quad_arc_vec_t *root, quad_arc_vec_t *A);
  /* Writes to {wr} a description of the map(s) that can be reached by
    {sym} and {onext} steps on the quad-edge structure from the root
    arcs {root->e[0..nr-1]}, where {nr = root->ne}.
    
    As a side effect, the procedure renumbers all edges in the reachable
    submap, with {quad_renumber_edges(root, A)}. */

void quad_read_map(FILE *rd, quad_arc_vec_t *root, quad_arc_vec_t *A);
  /* Reads from {rd} a description of a map in the format used by
    {write_map}, and builds the corresponding quad-edge structure in
    memory.
    
    The procedure allocates all the edge records needed to build the
    input map, and assigns them sequential numbers starting from 0. It
    also stores into the table {A} the base arc on each edge record,
    indexed by the edge number. It also reads from the file the root
    arcs of that map, and saves them in the table {root}. The tables
    {A} and {root} are (re)allocated or trimmed as needed, so that
    {A->ne} will be the number of edges in the map, and {root->ne}
    the number of roots. */

#endif
