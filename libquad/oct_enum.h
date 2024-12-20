#ifndef oct_enum_H
#define oct_enum_H

/* Element enumeration procedures for the general quad-edge structure. */
/* Last edited on 2024-12-05 10:39:34 by stolfi */

#define oct_enum_H_copyright \
  "Copyright © 1996, 2006 State University of Campinas (UNICAMP).\n\n" jslibs_copyright

#include <stdint.h>

#include <jslibs_copyright.h>
#include <bool.h>
#include <vec.h>

#include <oct.h>

/* ENUMERATING ARCS OF AN OCT-EDGE STRUCTURE
  
  For general concepts about data structure enumeration
  (especially /root item/, /step-func/, /visit-func/, /visit-list/),
  see {enum_orbits.h}. 
  
  Most enumeration procedures in this interface are specialized calls
  to the general {enum_cycle} and {enum_orbits} procedures from
  {enum_orbits.h}. However, the procedures below generally deal with
  {oct_arc_t} items, instead of amorphous {ref_t} items.
  
  Also, the visit-lists are arc vectors ({oct_arc_vec_t*}) instead of
  general pointer vectors ({ref_vec_t*}). */

typedef oct_arc_t oct_step_t(oct_arc_t a);
  /* A step-func for oct-edge structure enumeration, namely a function
    that maps arcs to arcs.
    
    Some of the procedures below, such as {enum_ring} and {enum_orbits},
    use arbitrary step-funcs provided by the client.  Other procedures
    use a fixed set of step-funcs, defined internally. */

typedef bool_t oct_visit_t(oct_arc_t p);
 /* A client-provided visit-function, that is called by an enumeration
  procedure to process each arc found in the enuemration.
 
  An {oct_visit_t} procedure must follow the general rules for
  visit-funcs explained in {enum_orbits.h}.

  In particular, a visit-func call that returns TRUE aborts the
  enumeration. The visit-func should not use long-jumps for that
  purpose. A {NULL} visit-func is treated as a no-op that always
  returns {FALSE}. */

/* GENERIC STRATIFIED ENUMERATION */

bool_t oct_enum
  ( oct_arc_vec_t root, 
    uint ns, 
    oct_step_t *step[], 
    oct_visit_t *visit[], 
    oct_arc_vec_t *vP
  );
  /* Equivalent to procedure {enum_items} in {enum_orbits.h} (q.v.),
    except that the items are {oct_arc_t}s instead of generic {ref_t}s.
    
    Namely, enumerates all arcs that can be reached from the arcs
    {root.e[0..root.ne-1]} by all possible sequences of calls to the
    step-funcs {step[0..ns-1]}. These are tried in order of increasing
    index, so that the procedure only attempts to use {step[i]} after
    exhausting all possibilities with steps {step[0..i-1]}.
    
    The procedure calls {visit[i](e)}, for {i} in {0..ns-1}, when it
    used {step[i]} to each a new arc {e} for the first time. It also
    calls {visit[ns](e)} when it gets to a root arc {e} that could not
    be reached from previous roots. Note that the list {step} must
    have {ns} elements, while {visit} must have {ns+1}. 
    
    A visit-proc may return TRUE to halt the enumeration with a TRUE
    result; otherwise the result is FALSE. A {visit[k]} that is NULL
    is equivalent to a visit-proc that does nothing and returns FALSE.
    
    If {vP} is not NULL, the procedure will store into {*vP} all
    visited elements. See the comments above and in {enum_orbits.h}
    for more details on the parameters {step}, {visit}, and {vP}. */

/* ENUMERATING ITEMS OF A CYCLE */

bool_t oct_enum_cycle
  ( oct_arc_t root, 
    oct_step_t *step, 
    oct_visit_t *visit, 
    oct_arc_vec_t *vP
  );
  /* Visits every arc that can be reached from the arc {root} by a
    sequence of zero or more applications of the {step} function.
    
    This is a special case of {oct_enum_orbits} with a single root and
    a single step-func. The {visit} function and the visit-list {vP},
    if not NULL, are used according to the general rules. */

/* ENUMERATING GENERAL ORBITS */

bool_t oct_enum_orbits
  ( oct_arc_vec_t root,
    uint ni, 
    oct_step_t *istep[], 
    uint no, 
    oct_step_t *ostep[], 
    oct_visit_t *visit,
    oct_arc_vec_t *vP
  );
  /* Enumerates all {istep}-orbits that can be reached fron the 
    {root} arcs by {ostep}-chains.
    
    More precisely, the procedure enumerates all arcs that can be
    reached from {root.e[0..root.ne-1]} by any sequence of steps in
    {istep[0..ni-1]} and/or {ostep[0..no-1]}. It then visits exactly
    one arc in each maximal subset that is connected by chains of
    {istep} functions alone. The {visit} function and the visit-list
    {vP}, if not NULL, are used according to the general rules. */

/* ORIENTATION-PRESERVING AND VIEW-PRESERVING ENUMERATIONS

  The following procedures will enumerate all arcs that can be reached
  from the arcs {root.e[0..root.ne-1]} by arbitrary combinations of
  *view-* and *orientation-preserving* walking functions (or,
  equivalently, by arbitrary sequences of {oct_sym} and {oct_onext}).
  
  A step-func {stp} is /view-preserving/ if it preserves the
  dual/primal character of the argument; that is, {a.M == a.stp.M} for
  every arc {a}, where {.M} is the map that {a} belongs to. In that
  case, every arc {a} found in the enumeration satisfies {a.M == r.M}
  for some root arc {r}. Thus, if all roots are on the same map, only
  arcs on that map will be visited. Of course, the arcs on the other
  map can be obtained by applying `{*}' to the visited arcs.
  
  A step-func {stp} is /orientation-preserving/ if its effect on any
  arc {a} is equivalent to a sequence of {oct_ostep} and/or
  {oct_fflip} moves, with an even number of the latter. If {Y} is an
  orientable connected component of the underlying manifold {X}, then
  any enumeration that starts from some root {r} on {Y} and uses only
  steps that are orientation-preserving within {Y} will only visit
  arcs that determine the same orientation on {Y} as {r}. In
  particular, if an arc {a} is visited, then {fflip(a)}, {vflip(a)},
  and {eflip(a)} are not; and vice-versa. In this case, the arcs with
  opposite orientation can be be obtained by applying any of those
  flip operations to every visited arc.

  On the other hand, if the the starting root {r} lies on a
  non-orientable connected component {Y} of {X}, then the same
  procedure may visit both {a} and one of {fflip(a)} and {vflip(a)},
  for some {a}. Ditto if the given root set contains two arcs that are
  oddly-connected to each other.

  Each procedure will visit exactly one arc {e} among all the
  enumerated arcs that are associated to the same map element of the
  designated type. */

bool_t oct_enum_arcs(oct_arc_vec_t root, oct_visit_t *visit, oct_arc_vec_t *vP);
  /* Visits all arcs reachable from {root} by view- and orientation-preserving
    paths.  Equivalent to {oct_enum} with {stp == (oct_onext,oct_sym)}
    and {vis == (visit,visit,visit)}. */

bool_t oct_enum_nodes(oct_arc_vec_t root, oct_visit_t *visit, oct_arc_vec_t *vP);
  /* Visit exactly one arc among all reachable arcs with the same
    origin node. Equivalent to {oct_enum_orbits} with {istep ==
    (oct_onext)} and {ostep == (oct_sym)}. */

bool_t oct_enum_edges(oct_arc_vec_t root, oct_visit_t *visit, oct_arc_vec_t *vP);
  /* Visit exactly one arc among all reachable arcs with the same
    undirected edge. Equivalent to {oct_enum_orbits} with {istep ==
    (oct_sym)} and {ostep == (oct_onext)}. */

bool_t oct_enum_faces(oct_arc_vec_t root, oct_visit_t *visit, oct_arc_vec_t *vP);
  /* Visit exactly one arc among all reachable arcs with the same left
    face. Equivalent to {oct_enum_orbits} with {istep == (oct_lnext)}
    and {ostep == (oct_sym)}. */
    
/* OCTET ENUMERATION */

bool_t oct_enum_octets(oct_arc_vec_t root, oct_visit_t *visit, oct_arc_vec_t *vP);
  /* Enumerates all arcs, of any view and orientation, that are
    reachable from the given {root} list. Then visits exactly one arc
    among each set of eight arcs which have the same dual pair of
    undirected edges.  Equivalent to {oct_enum_orbits} with {istep ==
    (oct_vflip,oct_fflip,oct_rot)} and {ostep == (oct_onext)}. */

#endif
