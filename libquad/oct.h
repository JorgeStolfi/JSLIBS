#ifndef oct_H
#define oct_H

/* The quad-edge data structure for orientable and non-orientable maps. */
/* Last edited on 2011-12-22 15:01:49 by stolfilocal */

#define oct_H_copyright \
  "Copyright © 1996, 2006 Institute of Computing, Unicamp."
  
/* THE (FULL) QUAD-EDGE STRUCTURE

  The quad-edge data structure encodes the topology of a 2-D map: 
  a graph drawn on a borderless surface of finite extent 
  (i.e. a compact 2-D manifold) in such a way that every face
  is a topological disk. For details, see

    [1] "Primitives for the Manipulation of General Subdivisions 
    and the Computation of Voronoi Diagrams".
    L. Guibas, J. Stolfi, ACM Transactions on Graphics (April 1985).

  The name {oct.h} is meant to distinguish this module, that
  implements the full data structure described in the article, from
  the numerous {quad.h} out there that only allow orientable maps and
  do not support the {Flip} operator.

  For readability, we often omit the {oct_} prefix from operations
  of this interface.

  AUTHORS

  This interface was originally created by J. Stolfi in April 1993. It
  was based on the orientable-map version {quad.h} implemented by Jim
  Roth (DEC CADM Advanced Group) in May 1986, which was the basis for
  the full Modula-3 interface {Oct.i3} by Rober M. Rosi in 1993. The
  latter was converted to C by J. Stolfi in 1996 and was substantially
  revised by him in January 2007. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <vec.h>
#include <bool.h>

/* MAPS

  A /map/, as represented by the quad-edge data strcuture,
  is a partition of some two-dimensional manifold {X} into 
  finitely many open {k}-dimensional balls, the /map elements/,
  satisfying certain incidence constraints.  The elements
  with dimension 0, 1, or 2 are called /vertices/, /edges/
  and /faces/ of the map. */

typedef struct oct_rep_t *oct_arc_t; 
  /* A value {e} of type {arc_t} specifies an /arc/ {e} of a quad-edge
    structure, which is an edge {E = e.edge} belonging to any of the
    two dual maps represented by the structure, taken with specific
    /longitudinal/ and /transversal orientations/.
    
    The longitudinal orientation of {e} is a choice among the two
    senses of motion along the edge {E}, as to which one is to be
    considered /forwards/ and which one is /backwards/. It allows us
    to distinguish the /destination/ vertex (reached by going forwards
    on {e}) from the /origin/ one (reached by going backwards).
    
    The transversal orientation of {e} is a choice among the two ways
    of crossing the edge {E}, as to which one is /leftwards/ and which
    one is /rightwards/. It allows us to distinguish the /left face/
    of {e} from its /right face/.
    
    The longitudinal and transversal orientations define a /circular
    orientation/ at any point {p} of {e}, allowing us to distinguish
    the /positive/ sense of turning around {p}
    (forwards-leftwards-backwards-rightwards) from the /negative/
    sense (forwards-rightwards-backwards-leftwards). */

/* ARC ORIENTATION OPERATORS

  The operators in this section change the orientations and primal/dual
  character of an arc {e}, preserving the primal and dual edge pair. */

oct_arc_t oct_nop(oct_arc_t e);
  /* The identity or `no operation' function that returns {e}
    itself; it is provided here for completeness. */ 

oct_arc_t oct_sym(oct_arc_t e);
  /* The {sym} operator reverses the longitudinal and transversal
    orientations. Intuitively, it rotates {e} by a half-turn. Thus
    {sym(e) != e} but {sym(sym(e)) = e}, for any {e}. */
     
oct_arc_t oct_fflip(oct_arc_t e);
  /* Reverses the transversal orientation of {e}, preserving its
    longitudinal orientation and underlying edge. The result is always
    distinct from {e}, {rot(e)}, {sym(e)}, and {tor(e)}; but
    {fflip(fflip(e)) == e} and {fflip(rot(e)) == tor(fflip(e))}. Thus,
    {fflip} also reverses the local circular orientation of {e}. */

oct_arc_t oct_vflip(oct_arc_t e);
  /* Reverses the longitudinal orientation of {e}, preserving its
    transversal orientation and primal/dual character. Thus
    {vflip(fflip(e)) = sym(e) = fflip(vflip(e))}. It
    also reverses the local circular orientation, so
    {vflip(rot(e)) == tor(vflip(e))}. */

/* DUALITY OPERATORS

  Actually, at any moment the quad-edge data structure represents
  *two* maps on {X}, each the dual of the other.
  
  The two maps are disjoints, and there is a one-to-one /duality/
  correspondence {*} between their elements such that, for any element
  {A} of either map, the dimensions of {A} and {A*} add to 2.
  Moreover, {A*} intersects {A} at precisely one point, and is
  disjoint from {B} or {B*}, for any other element {B!= A} with the
  same dimension as {A}. */

oct_arc_t oct_rot(oct_arc_t e);
oct_arc_t oct_tor(oct_arc_t e);
  /* The {rot} operation returns an arc {f} on the dual edge {F} of
    the edge {E} underlying {e}, oriented so that {f} crosses {e}
    leftwards and {e} crosses {f} rightwards. Intuitively, {rot}
    rotates {e} by a positive quarter-turn. Thus, {rot(e) != e},
    {rot(e) != sym(e)}, but {rot(rot(e)) == sym(e)}, for any {e}.
    
    The {tor} operator is the inverse of {rot}: it rotates {e} by a
    negative quarter-turn. It is equivalent to {rot(sym(e))}.
    
    Both {rot} and {tor} preserve the circular orientation 
    of the arc itself. */
  
oct_arc_t oct_dual(oct_arc_t e);
  /*  The /left-dual/ of the arc {e}: the arc {f} on the dual edge
    {e.edge*}, oriented so that so that {f} crosses {e} leftwards, and
    {e} crosses {f} leftwards.
    
    Said another way, the longitudinal orientation of {dual(e)} agrees
    with the transversal orientation of {e}, and vice-versa. It is
    equivalent to {fflip(rot(e))}, {tor(fflip(e))}, {vflip(tor(e))}, and
    {rot(vflip(e))}.
    
    Note that {e} and {dual(e)} define opposite circular orientations
    on any disc neighborhood that contains both. Note also that
    {dual(dual(e)) == e} for any {e}. */ 

oct_arc_t oct_duar(oct_arc_t e);
  /* The /right-dual/ of {e}: the arc {f} on the dual edge {e.edge*},
    oriented so that {f} crosses {e} rightwards, and {e} crosses {f}
    rightwards.
    
    Said another way, the longitudinal orientation of {duar(e)}
    *disagrees* with the transversal orientation of {e}, and
    vice-versa. It is equivalent to {fflip(tor(e))}, {rot(fflip(e))},
    {vflip(rot(e))}, and {tor(vflip(e))}.
    
    Note that {e} and {duar(e)} define opposite circular orientations on
    their common point. Note also that {duar(duar(e)) == e} for any {e}. */

/* TUMBLING OPERATORS 
 
  The /tumbling operators/ are the eight functions {nop}, {rot},
  {sym}, {tor}, {fflip}, {vflip}, {dual}, and {duar}. This set is
  closed under composition: any sequence of tumblings is equivalent to
  a single tumbling. Namely, the composition {G(F(e)) = e.F.G} is
  given by the table

            | G                                           
            | nop   sym    fflip vflip  rot   tor    dual  duar
    --------|---------------------------------------------------
    F nop   | nop   sym    fflip vflip  rot   tor    dual  duar
      sym   | sym   nop    vflip fflip  tor   rot    duar  dual
            |                                 
      fflip | fflip vflip  nop   sym    DUAR  DUAL   TOR   ROT
      vflip | vflip fflip  sym   nop    DUAL  DUAR   ROT   TOR
            |                                 
      rot   | rot   tor    DUAL  DUAR   sym   nop    VFLIP FFLIP 
      tor   | tor   rot    DUAR  DUAL   nop   sym    FFLIP VFLIP
            |                                 
      dual  | dual  duar   ROT   TOR    FFLIP VFLIP  nop   sym
      duar  | duar  dual   TOR   ROT    VFLIP FFLIP  sym   nop  

  Note that the operators {sym,rot,tor} commute, {sym,fflip,vflip}
  commute, and {sym,dual,duar} commute; but {fflip,vflip} do not
  commute with {rot,tor}, and neither pair commutes with {dual,duar}.
  
  Note also that the tumblings under composition constitute a group,
  isomorphic to the symmetry group of a non-cubical square prism. */

oct_arc_t oct_rot_fflip(oct_arc_t e, int r, int f);
  /* This generic operator changes the orientation of {e}
    by applying {rot} to {e}, {r} times, then applying
    {fflip} to the result, {f} times.  
    
    For negative values of {r} and/or {f}, the result is defined by
    the identities {rot_fflip(e,-r,f) = rot_fflip(e,3*r,f)} and {rot_fflip(e,r,-f) =
    rot_fflip(e,r,f)}. Thus, for example,
    
      {rot_fflip(e,1,0) == rot(e)}, 
      {rot_fflip(e,2,0) == sym(e)},
      {rot_fflip(e,3,0) == tor(e)},
      {rot_fflip(e,0,1) == fflip(e)}, 
      {rot_fflip(e,-5,7) = fflip(tor(e))},
      {rot_fflip(e,4,-6) = e},
      
    and so on. */

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
  or right ('{r}') bank of the arc. The result of the operator is
  another arc whose circular orientation also agrees with {C}, by the
  same rule.

  Note that the same edge may appear twice in the list {S}. In such
  cases, the two occurences will have opposite longitudinal
  orientations (if the pivot is a vertex) or opposite transversal
  orientation (if the pivot is a face). In either case, the other
  orientation may or may not be reversed. On the other hand, in
  ordinary quad-edge structures the list {S} of {e} will never include
  any of its four dual versions ({rot(e)} and its
  {fflip/vflip} variants). */

oct_arc_t oct_onext(oct_arc_t e);
oct_arc_t oct_oprev(oct_arc_t e);
 /* The next and previous arc, respectively,
   among the arcs with same origin as {e}.  */

oct_arc_t oct_dnext(oct_arc_t e);
oct_arc_t oct_dprev(oct_arc_t e);
 /* The next and previous arc, respectively,
   among the arcs with same destination as {e}.  */

oct_arc_t oct_lnext(oct_arc_t e);
oct_arc_t oct_lprev(oct_arc_t e);
 /* The next and previous arc, respectively,
   among the arcs with same left face as {e}.  */

oct_arc_t oct_rnext(oct_arc_t e);
oct_arc_t oct_rprev(oct_arc_t e);
 /* The next and previous arc, respectively,
   among the arcs with same right face as {e}.  */

oct_arc_t oct_walk(oct_arc_t e, int r, int n);
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
    
/* !!! Should we define also {eflip(a) = fflip(onext(a))}? !!! */

/* COUNTING INCIDENCES */

uint oct_odegree(oct_arc_t e);
  /* The number of arcs that leave the origin vertex of arc {e}. */
  
uint oct_ddegree(oct_arc_t e);
  /* The number of arcs that arrive at the destination vertex of arc {e}. */
  
uint oct_ldegree(oct_arc_t e);
  /* The number of arcs that have the same left face as {e}. */

uint oct_rdegree(oct_arc_t e);
  /* The number of arcs that have the same right face as {e}. */

/* TOPOLOGY-CHANGING OPERATORS */

oct_arc_t oct_make_edge(void);
  /* Adds to the quad-edge structure a new connected component,
    consisting of a map on the sphere with a single edge, a 
    single vertex, and two faces.  The edge number is set to zero;
    use {oct_set_edge_num} to change it.  */

void oct_destroy_edge(oct_arc_t e);
  /* The arc {e} must belong to a connected component of the map with
    a single edge, in any topology. Reclaims the storage used by that
    component. After this call, passing {e} or any of its turned
    variants to any procedure in this interface will have
    unpredictable results. */

void oct_splice(oct_arc_t a, oct_arc_t b);
  /* Performs a splicing operation on the map, by exchanging 
    {onext(a)} and {onext(b)}, as well as {lprev(a)} and {lprev(b)}.
    
    This procedure should not be called if {b == fflip(onext(a))},
    or if {b} can be reached from {a} by a sequence of {rot}
    and {onext} operations with an odd number of {rot}s. */

/* EDGE RECORDS */

typedef struct oct_edge_rec_t *oct_edge_t;
  /* An edge record, representing an undirected and unoriented edge {E}
    of the primal map together with the corresponding edge {F} of
    the dual map. It is shared by the eight arcs that consist of {E} or {F}
    taken with all possible logitudinal and transversal orientations. */
  
oct_edge_t oct_edge(oct_arc_t e);
  /* Obtains the edge record of an arc reference {e}. Satisfies
    {edge(rot(e)) == edge(fflip(e)) == edge(e)}.  */

/* ARC AND EDGE NUMBERS */

uint64_t oct_edge_num(oct_edge_t E);
void oct_set_edge_num(oct_edge_t E, uint64_t num);
  /* These functions get and set the current number of an edge {E}.
    The number must be limited to 61 bits (i.e. in {0..2^61-1}). */

uint64_t oct_arc_num(oct_arc_t *e);
  /* Returns a numeric identifier for the arc {e}, consisting of
    {8*edge_num(edge(e)) + tumble_code(e)}. Thus, if edge numbers are
    unique in a given set of edges, then the arc numbers will be
    unique among the set of all arcs on those edges.
    
    Other sorts of arc numbers can be obtained by combining the edge
    number and the orientation bits, in a similar way. For instance,
    the formula {2*edge_num(edge(e)) + lon_bit(e)} gives unique
    numbers to all longitudinally directed primal edges. */

/* TUMBLE CODE

  A value {e} of type {arc_t} consists of a pointer
  {p = edge(e)} to an /edge record/, of type {edge_rec_t}, 
  and a three-bit /tumble code/ {t = tumble_code(e)}.
  
  The pointer {p} and the tumble code {t} are packed together as a single
  {void *} value. This trick assumes that the address of any edge
  record is a multiple of 8 bytes (64 bits), so {t} can be stored in
  its lowest three bits, without ambiguity. Thus, although an
  {arc_t} is formally an address, it should not be used as such,
  since it generally points to some random byte inside the edge
  record. */

typedef unsigned char oct_bits_t;
   /* A data type used to hold a few bits (usually at the low-order end). */
    
oct_bits_t oct_tumble_code(oct_arc_t e);
  /* Returns a three-bit code (an integer in {0..7}) that identifies
    the arc {e} among the eight arcs with the same edge record. Two
    arcs {e,f} are identical if and only if {edge(e)==edge(f)} and
    {tumble_code(e)==tumble_code(f)}. */
  
oct_arc_t oct_orient(oct_edge_t E, oct_bits_t t);
  /* Returns the unique arc that has edge record {E} and tumble code {t}. 
    Only the three lower-order bits of {t} are used.
    
    This procedure can be convenient in special situations, like
    reading or writing a quad-edge structure, or enumerating all arcs
    with the same edge. It is not very useful in general because the
    meaning of the tumble codes is implementation-dependent. In
    particular, {orient(edge(e),tumble_code(f))} is hardly sensible
    if {edge(f)!=edge(e)}. */

/* ORIENTATION BITS
  
  The tumble code of an arc {e} can be mapped to three Boolean
  attributes: the /longitudinal orientation bit/ {lon_bit(e)}, the
  /transversal orientation bit/ {trn_bit(e)}, and the /dual-primal
  bit/ {dop_bit(e)}. These three bits completely determine the tumble
  code (in some implementation-dependent way). */
  
oct_bits_t oct_lon_bit(oct_arc_t e);
  /*  The /longitudinal orientation bit/ of {e}. This attribute is
    reversed by {sym} and {vflip}, but preserved by {fflip}. It can be
    used to distinguish the pair {e,fflip(e)} from {sym(e),vflip(e)}. */
  
oct_bits_t oct_trn_bit(oct_arc_t e); 
  /* The /transversal orientation bit/ of {e}. This attribute is
    reversed by {sym} and {fflip}, but preserved by {vflip}. It can be
    used to distinguish the pair {e,vflip(e)} from the pair
    {sym(e),fflip(e)}. */
  
oct_bits_t oct_dop_bit(oct_arc_t e); 
  /* The /dual-primal bit/ of {e}. This attribute is reversed by {rot}, {tor}, {dual}
    and {duar}, but preserved by {fflip}, {vflip} and {sym}. Therefore, it
    can distinguish between the arcs of the edge {E} ({e,sym(e),fflip(e),vflip(e)})
    and those of its dual edge ({rot(e),tor(e),dual(e),duar(e)}). */

oct_bits_t oct_cir_bit(oct_arc_t e);
  /* The /circular orientation bit/ {cir_bit(e)} is defined as the
    exclusive-or of {lon_bit(e)} and {trn_bit(e)}. This attribute
    is preserved by {rot,tor,sym}, but reversed by {fflip,vflip,dual,duar}.
    Threfore, it can distinguish between the four arcs {e,rot(e),sym(e),tor(e)},
    that imply the the same circular orientation around a point, 
    from {fflip(e),vflip(e),dual(e),duar(e)}, that imply the opposite sense. */

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
  
  There is also no constraint on the {lon_bit}s of {e} and {rot(e)}.
  In fact, if {lon_bit(e)==lon_bit(rot(e))} for some arc {e}, then one
  must have {lon_bit(f)!=lon_bit(rot(f))} for the arc {f==rot(e)}.
  
  If the surface is orientable, the quad-edge structure can be set up
  so that the {cir_bit}s of {e}, {e.onext}, and {e.lnext} are
  identical for all arcs {e}. In that case, the {fflip} and {vflip}
  operators are trivial, and one could use the quad-edge structure
  instead. (However, the {fflip} operator may be useful even when
  working with orientable maps; for example, it makes it possible to
  produce two identical copies of an orientable but asymmetric map,
  and then join one copy with a mirror image of the other, in {O(1)}
  operations.) */

/* NULL ARCS */

#define oct_arc_NULL ((oct_arc_t)NULL)
  /* A special {arc_t} value that means `no such arc'. */

bool_t oct_is_null(oct_arc_t e);
  /* TRUE if {e} is a null {arc_t}, equivalent to {oct_edge(e)==NULL}.
    Note that this is a stronger test than {e == NULL} or {e ==
    oct_arc_NULL} because one can have tumbled versions of the {NULL}
    edge. */

/* VECTORS OF ARCS AND EDGES */

vec_typedef(oct_arc_vec_t,oct_arc_vec,oct_arc_t);
  /* An {oct_arc_vec_t} is a vector of {oct_arc_t}s. */

vec_typedef(oct_edge_vec_t,oct_edge_vec,oct_edge_t);
  /* An {oct_edge_vec_t} is a vector of {oct_edge_t}s. */
  
/* ARC INPUT/OUTPUT  */

void oct_write_arc(FILE *wr, oct_arc_t e, int width);
  /* Writes the arc {e} to {wr} in the format "{num}:{t}" where {num}
    is {edge_num(edge(e))}, and {t} is {tumble_code(e)} in binary.
    The {num} is padded with spaces to have at least {width}
    characters. */

oct_arc_t oct_read_arc(FILE *rd, oct_arc_vec_t *et);
  /* Parses from file {rd} a description of an arc, in the format
    "{num}:{t}" where {num} is an edge number and {t} is a tumble
    code. Requires an edge table {et}, which is a vector of arcs such
    that {et->e[i]} is any arc on an edge edge {E} with {E->num == i}.
    The resulting arc {e} has {oct_edge(e) == et->e[num]} and
    {oct_tumblecode(e) == t}. */

/* MAP INPUT/OUTPUT  */

void oct_write_map(FILE *wr, oct_arc_vec_t *root, oct_arc_vec_t *et);
  /* Writes to {wr} a description of the map(s) that can be reached by
    {sym} and {onext} steps on the quad-edge structure from the root
    arcs {root->e[0..nr-1]}, where {nr = root->ne}.
    
    As a side effect, the procedure renumbers all octets (primal/dual
    edge pairs) in the reachable submap, sequentially from 0.
    If {et} is not NULL, stores into {et->e[0..et->ne-1]} one arc
    out of every reachable octet, indexed by the edge number.
    The table {*et} is trimmed or expanded by the procedure
    as needed. */

void oct_read_map(FILE *rd, oct_arc_vec_t *root, oct_arc_vec_t *et);
  /* Reads from {rd} a description of a map in the format used by
    {write_map}, and builds the corresponding quad-edge structure in
    memory.
    
    The procedure allocates all the edge records needed to build the
    input map, and assigns them sequential numbers starting from 0. It
    also stores into the table {et} the base arc on each edge record,
    indexed by the edge number. It also reads from the file the root
    arcs of that map, and saves them in the table {root}. The tables
    {et} and {root} are (re)allocated or trimmed as needed, so that
    {et->ne} will be the number of edges in the map, and {root->ne}
    the number of roots. */

#endif
