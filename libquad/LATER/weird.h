#ifndef weird_H
#define weird_H

/* A data structure inspired on Weiler's radial edge for 2D pseudomanifolds. */
/* Last edited on 2024-07-01 18:55:32 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>

#define weird_H_copyright \
  "Copyright © 2024 State University of Campinas (UNICAMP).\n\n" jslibs_copyright

#include <jslibs_copyright.h>
#include <bool.h>
#include <vec.h>

/* THE WEIRD STRUCTURE 

  The /weird/ data structure, inspired on Weiler's "radia edge"
  structure, is meant to encode the topology and geometry of a /mesh/ consisting of a
  finite number of /elements/ which are open polygons (the /faces/) open
  straight line segments (the /edges/), and points (the /vertices/).
  with certain incidence and adjacency conditions.
  
  TOPOLOGY: FLAGS AND FLAG OPERATORS
  
  The topology part of a weird structure {W} is a quintuple {T(W) =
  (F,.sym,.enext,.fnext)} where {F = F(W)} is a finite set -- the
  /flags/ of {W} -- and {.sym}, {.enext}, and {.fnext} --
  the /operators/ of {W} -- are permutations of {F}, with properties
  specified below.
  
  The {.enext} and {.fnext} operators may be arbitrary permutations of
  the flags. Their inverses will be denoted {.eprev}, and {.fprev},
  respectively. The {.sym} operator must be a proper involution (meaning
  that {f.sym != f} and {f.sym.sym = f} for any flag {f}), and it must
  satisfy the properties {.enext.sym} = {.sym.eprev} and {.fnext.sym =
  .sym.fprev}. It folows that {sym.enext = .eprev.sym} and {.sym.fnext =
  .fprev.sym} and {.sym.vnext = .vprev.sym}.
  
  (The topology can be seen as derived from a /3-G-map/ as defined by
  Lienhardt [1a] and Guedes [1b], where {.sym} is the operator
  {\phi_0\phi_3}, {.fnext} is {\phi_0\phi_1}, and {.enext} is
  {\phi_2\phi_3}. This in turn is a special case of the GEM
  (graph-encoded manifold) structure described and studied by Lins [2];
  which is a special case of the GEM structure defined by Montagner and
  Stolfi [3]. See also to the facet-edge structure of Laslo and Dobkin
  [4].
  
  More specifically, the topology of a weird structure is derived
  from an *orientable* 3-G-map, whose operator graph is
  bipartite. In this case every G-map operator swiches side of the
  bipartition, whereas the composite operators of the weird structure
  stay on the same side. Then the flags {F} of the weird structure are
  just one side of the bipartition.)
  
  (Just as we can go from the quad-edge structure to the oct-edge structure,
  we could replace {F} by two copies {F+,F-}, with {.sym}, {.enext}, and {.fnext}
  and add a fifth operator {.flip} to the weird structure with the
  requirement that {.op.flip = .flip.op^{-1}} for every operator {op} among the 
  other three.  How it relates to the other four is to be defined.)
  
  RINGS AND STARS
    
  Each topological operator {.op} of theweird structure  determines a partition of the flags into
  one or more orbits, each orbit consisting of {f.op^k} for every
  integer {k}.
    
  The orbits of {.sym} are pairs of flags.  Each pair is called a /topological edge/
    
  Each orbit of {.fnext} 
  
  
  But
  {f.cflip.enext^k = f.eprev^k.cflip}; because of W7, the orbits of {f}
  and {f.cflip} are disjoint.
  
  Similarly, let's define {.fnext = .vflip.eflip}; its inverse {.fprev} will be {.eflip.vflip}.
  
  Also let's define {.vnext = .fflip.eflip}; its inverse {.vprev} will be {.eflip.fflip}.
  
  A weird structure is reqired to satisfy two more topological axioms:
    
    W8. For any {f} and any {k}, {f.enext^k = (f(.vflip.eflip)^k} is distinct from {f.cflip}.
  
  An orbit of 
  
  TOPOLOGICAL INTERPRETATION
  
  A weird data structure {W} can be interpreted as a topological space
  {T(S}} as follows:
  
    For each quad {q} of {W}, choose one flag {f(q)} to be its /base flag/
    Any flag {g} in {q} can be specified by a pair of bits {o,r} such that 
    {g = f(q).vflip^o.cflip^r = f(q).cflip^r.vflip^o}. This bit pair is the /orientation/
    of the flag {g} (with respect to the chosen base flag of its quad).  Then we can denote 
    flag {g} by {q[o,r]}.
    
    For each quad {q} of {W}, define {T(q)} to be a copy of the triangle
    {(0,0),(0,1),(1,0)} of {\RR^2}.
    
    For any quad {q} and any orientation {o,r}, let {q'[o',r']} be the 
    flag {q[o,r].eflip}. 
    
  GEOMETRY: 
  
  A flag
  can  be imagined as 
  arcs of the same edge are said to be /opposite/ of each other. The
  (/oriented/) /border/ of each face is a circular list of alternating
  arcs and vertices, where each flag is oriented in the direction that
  matches counterclockwise turning on the face. Each flag may appear only
  once in this list, and only on one face (but an flag and its opposite may
  both appear on the same face).
  
  The half-edge structure consists of /flag records/, each representing,
  by definition, one flag of the mesh. The structure defines two
  fundamental /walking operators/ on ac records, {.cflip} and {.lnext}.
  The {.cflip} operator is an involution of set of all flag records, with
  no fixed points: {a.cflip != a} and {a.cflip.cflip = a} for every flag {a}. The
  records {a} and {a.cflip} are interpreted as the oppositely oriented
  arcs on the same edge of the mesh.
  
  The {.lnext} operator is an arbirary permutation of the flag records,
  meaning that {a.lnext != b.lnext} if {a != b}. Each cycle of this
  permutation is interpreted as the list of arcs on the border of a face
  of the mesh.
  
  The {.oprev} operator is defined as the composition {.cflip.lnext},
  applied in that order. Each cycle of {.oprev} corresponds to a unique
  vertex of the mesh, and consists of the arcs that leave that vertex,
  in /clockwise/ order around the vertex. Again, an flag {a} may appear
  only once in only one of these cycles; but {a} and {a.cflip} may both
  appear on the same cycle, signifying that the corresponding edge is a
  loop with both endpoints on the same vertex.
  
  Several other walking operators can be defined in terms of these, such as
  {.lprev}, {.onext}, {.rnext}, {.rprev}, {.dprev}, {.dnext}. See below.
  
  Each flag record {a} also contains a non-negative integer /flag id/
  {a.aid} that can be chosen by the client of this interface. The flag
  ids of {a} and {a.cflip} are always {2*a.eid} and {2*a.eid+1} or
  vice-versa, where {a.eid=a.cflip.eid} is a client-chosen id number for
  the underlying undirected edge. The integers {a.aid} and {a.eid} can
  be used as indices in external tables to attach arbitrary properties
  to each flag or edhe of the mesh.
  
  The half-edge data structure proper has no explicit records for faces
  and vertices. Clients of this interface who need to attach properties
  to those elements of the mesh (such as vertex coordinates or face
  geometry) may define their own data records for them. In that case,
  they should probably set up tables {left[a.aid]} that gives the record
  for the face that has {a} on its border, and/or {org[a.aid]} for the
  record of the vertex that is the origin of flag {a}. See {weird_enum_faces}
  and {weird_enum_vertices} below. */

typedef struct weird_rep_t *weird_flag_t;
  /* An /flag reference/, a pointer to a record of the half-edge 
    structure representing an /flag/ (directed edge) of the mesh. */

vec_typedef(weird_flag_vec_t, weird_flag_vec, weird_flag_t);
  /* An extensible vector of flag references. */

/* FAST WALKING OPERATORS */

weird_flag_t weird_ref(weird_flag_t a);
  /* Same edge in opposite orientation.  Takes constant time. */
  
weird_flag_t weird_lnext(weird_flag_t a);
  /* Next flag /counterclockwise/ with the same left face. Takes constant time. */
 
weird_flag_t weird_rnext(weird_flag_t a);
  /* The next flag /counterclockwise/ with the same right face,
    namely {a.cflip.lnext.cflip}. Takes constant time. */
       
weird_flag_t weird_oprev(weird_flag_t a);
  /* The next flag /clockwise/ with the same origin vertex,
    namely {a.cflip.lnext}. Takes constant time. */
  
weird_flag_t weird_dprev(weird_flag_t a);
  /* he next flag /clockwise/ with the same destination vertex
    namely {a.lnext.cflip}. Takes constant time. */

/* SLOW WALKING OPERATORS */

weird_flag_t weird_lprev(weird_flag_t a);
  /* Inverse of {.lnext}: the next flag /clockwise/ with the same left face.
    Takes time proportional to the degree of that face. */

weird_flag_t weird_rprev(weird_flag_t a);
  /* Inverse of {.rnext}: next flag /clockwise/ with same right face.
    Takes time proportional to the degree of that face. */
      
weird_flag_t weird_onext(weird_flag_t a);
  /* Inverse of {.oprev}: the next flag /counterclockwise/ with the same origin vertex.
    Takes time proportional to the degree of that vertex. */
   
weird_flag_t weird_dnext(weird_flag_t a);
  /* Inverse of {.dprev}: the next flag /counterclockwise/ with the same destination vertex.
    Takes time proportional to the degree of that vertex. */

/* FLAG DIRECTION BIT 

  The /direction bit/ of an flag, either 0 or 1, distinguishes its from
  its opposite. Thus {a} and {a.cflip} have direction bits 0 and 1 or
  vice-versa. An flag with direction bit 0 is said to be the /base flag/
  of the underlying edge. */

typedef uint8_t weird_dir_bit_t;
   /* A data type used to hold a few bits (usually at the low-order end). */

weird_dir_bit_t weird_dir_bit(weird_flag_t a);
  /* The direction bit of flag {a}. */

weird_flag_t weird_base_flag(weird_flag_t a);
  /* The base flag with same edge as {a}, with direction bit 0. */

/* UNORIENTED EDGES */

typedef struct weird_edge_rec_t *weird_edge_t;
  /* Reference to an undirected and unoriented edge {ed}
    of the mesh.  It is shared by the two arcs that consist of {ed}
    taken with both directions. */

vec_typedef(weird_edge_vec_t, weird_edge_vec, weird_edge_t);
  /* An externsible vector od {half_edge_t} references. */

weird_edge_t weird_edge(weird_flag_t a);
  /* Obtains the edge reference of an flag reference {a}. Satisfies
    {a.cflip.edge == a.edge}.  */

weird_flag_t weird_orient(weird_edge_t ed, weird_dir_bit_t db);
  /* Returns the flag {a} that has {weird_edge(a) = ed} and {had_dir_bit(a) = db}.  */

/* EDGE IDENTIFIERS 
  
  The structure stores a identifying non-negative /id number/ for each
  undirected edge of the mesh. See {weird_set_edge_id} and
  {weird_renumber_edges} below for how that number gets defined.
  
  The edge identifier is often meant to identify each edge uniquely,
  e.g. for debugging or table indexing; but this interface does not 
  guarantee that by itself. */

typedef uint64_t weird_edge_id_t;  /* An edge ID number. */

#define weird_edge_id_MAX (weird_flag_id_MAX/2)
  /* Max edge identifier (so that the flag identifier does not overflow). */

weird_edge_id_t weird_edge_id(weird_flag_t a);
  /* Returns the identifier of the undirected edge underlying flag {a},
    namely {aid/2} where {aid} is the flag's id. */

typedef uint64_t weird_edge_count_t;  /* A count of edges in a data structure, table, etc.. */

#define weird_edge_count_MAX (weird_flag_count_MAX/2)
  /* Max number of edges in a structure. */

/* FLAG IDENTIFIERS 

  The /flag identifier/ of an flag {a} is {2*eid + db} where
  {eid} is the identifying number of the undirected edge, 
  and {db} is the direction bit. */

typedef uint64_t weird_flag_id_t;   /* An flag ID number. */

weird_flag_id_t weird_flag_id(weird_flag_t a);
  /* Returns the identifier of flag {a}. */

#define weird_flag_id_MAX (UINT64_MAX)
  /* Max flag identifier. */

typedef uint64_t weird_flag_count_t;   /* A count of arcs in a data structure, table, etc.. */

#define weird_flag_count_MAX (((uint64_t)1024)*1024*1024)
  /* Max number of arcs in a structure for sane table allocation. */

/* BUILDING AND MODIFICATION */
  
weird_flag_t weird_make_stick(weird_edge_id_t eid);
  /* Creates a new half-edge structure consisting of a pair of
    {weird_flag_t} records, {a,b}, representing a mesh with sperical
    topology, one non-loop edge, one face, and two distinct vertices.
    The records will have {a.cflip=b}, {b.cflip=a}, {a.lnext=b},
    {b.lnext=a}. The flag ids will be {2*eid} and {2*eid+1}. */

weird_flag_t weird_make_loop(weird_edge_id_t eid);
  /* Creates a new half-edge structure consisting of a pair of
    {weird_flag_t} records, {a,b}, representing a mesh with spherical
    topology, one loop edge, one vertex, and two distinct faces.
    The records will have {a.cflip=b}, {b.cflip=a}, {a.lnext=a},
    {b.lnext=b}. The flag ids will be {2*eid} and {2*eid+1}. */

void weird_splice(weird_flag_t a, weird_flag_t b);
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
 
void weird_set_lnext(weird_flag_t a, weird_flag_t b);
  /* Sets {a.lnext} to be {b}.  This operation may break the topological consistency.
    The caller must make sure that eventually the consistency is restored, namely that 
    the {.lnext} operator is a permutation of all the {weird_flag_t} records. */

void weird_set_edge_id(weird_flag_t a, weird_edge_id_t eid);
  /* Sets the identifiers of {a} and {a.cflip} to {2*eid} and {2*eid+1}. */
  
/* DEBUGGING */

void weird_check_topology(weird_edge_count_t ne, weird_flag_t a[], weird_edge_id_t eid0, bool_t verbose);
  /* Expects that the vector {a[0..ne-1]} has one {weird_flag_t} record
    out of every edge of the mesh.  
    Checks that the ids of flag {a[ke]} and {a[ke].cflip} are {2*(eid0+ke)} and {2(eid0+ke)+1}
    for all {ke} in {0..ne-1}.  Checks that the {.lnext} operator is a permutation of 
    those records and their opposites, and that the {.cflip} operator is an involution 
    on those records and their syms, with no fixed point. */

#endif
