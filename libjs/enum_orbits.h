#ifndef enum_orbits_H
#define enum_orbits_H

/* General stratified enumeration routines for groups and directed graphs. */
/* Last edited on 2016-12-07 07:20:45 by jstolfi */

#include <stdint.h>

#include <vec.h>
#include <bool.h>

/* DIRECTED GRAPH ENUMERATION 

  For this interface, an /enumeration procedure/ locates all items
  that can be reached from a given set of /root items/ by a chain of
  zero or more /steps/, in all possible combinations; and /visits/
  each item found in this process.

  Visiting an item means passing it to the caller for arbitrary
  processing (see {VISIT-FUNCTIONS} below) and/or saving it into an
  array provided by the caller (see {VISITATION LISTS} below). An
  enumeration procedure never visits the same item more than once.

  The items here are represented by {void *} pointers ({ref_t}
  values). The enumeration procedures never de-reference those
  pointers, so they need not be valid addresses.

  The enumeration procedures below do not modify the data structure or
  any global storage (except possibly through the the client's
  visit-func). Therefore, concurrent calls (e.g. in multi-threaded
  programming or enumeration-within-enumeration) do not
  directly interfere with each other. */

/* STEP-FUNCTIONS

  Each step of the enumeration consists of following some
  /step-func/, a function that maps items to items, provided by the
  caller. A step-func needs not be one-to-one. */

typedef ref_t enum_step_t(ref_t p);
 /* A step-func, namely a function from items to items. */

/* VISIT-FUNCTIONS

  An enumeration procedure takes as parameter the address of a
  /visit-func/ (of type {enum_visit_t}), and applies that function
  to each item found in the enumeration.

  The visit-func may return TRUE to abort the enumeration, or FALSE
  to allow it to continue. The enumeration procedure itself will
  return TRUE iff any visit-func call returned TRUE.

  A visit-func should not try to abort the enumeration with a long
  jump, or by killing its thread, since doing so will not give
  the enumeration procedure a chance to free any
  working storage that it may have allocated.  .

  A {NULL} visit-func is treated as a no-op that always returns
  {FALSE}. */

typedef bool_t enum_visit_t(ref_t p);
 /* A client-provided visitation procedure. */

/* VISIT-LISTS

  An enumeration procedure also takes a parameter {vP} of type
  {ref_vec_t*}. If {vP} is not {NULL}, the procedure will store into
  {vP} each item that it visits, in order of visitation.

  Typically, a client that uses this feature will set {vP} to the
  address of an empty {ref_vec_t}, with {vP->e == NULL} and {vP->ne
  == 0}. (Note that this is different from passing {vP == NULL}.) In
  any case, the new entries will be inserted *after* any entries
  that are already in {*vP}, namely starting at position {vP->e[vP->ne]}.

  Visited items are stored in {vP} even if the visit-func is NULL, or
  if the enumeration was aborted because a call to the visit-func
  returned TRUE.

  Upon exit, the count {vP->ne} will have been incremented by the
  number of items visited. The storage area {*(vP->e)} is expanded
  as needed, and trimmed at the end of the enumeration to hold
  precisely {vP->ne} items. */

/* SINGLE-CYCLE ENUMERATION */

bool_t enum_cycle(ref_t p, enum_step_t *step, enum_visit_t *visit, ref_vec_t *vP);
  /* Enumerates every item {q} reachable from {p} by zero or more
    calls to the step function {step}, and calls the visit-function
    {visit} on each, starting with {p} itself. See the comments
    above for details on {step}, {visit}, and {vP}.
    
    The function {step} must eventually lead back to {p}, otherwise
    the procedure will loop forever. */

/* ORBIT ENUMERATION */

bool_t enum_orbits
  ( ref_vec_t root, 
    uint32_t ni, 
    enum_step_t *istep[], 
    uint32_t no, 
    enum_step_t *ostep[], 
    enum_visit_t *visit,
    ref_vec_t *vP
  );
  /* Enumerates all {istep}-orbits that can be reached by {ostep}-chains
    from {root} items.
    
    More precisely, enumerates all items that can be reached from the
    items {root.e[0..root.ne-1]} by any sequence of steps in
    {istep[0..ni-1]} and/or {ostep[0..no-1]}. Then visits exactly one
    item in each maximal subset that is connected by chains of
    {istep[0..ni-1]} alone.
    
    The description above is accurate only if the step-funcs
    {istep[0..ni-1]} are one-to-one. A more general description is:
    whenever the procedure reaches a new item, it will exhaust all
    {istep} paths from that item, before trying any {ostep} function
    or looking at the next root.
    
    The visited elements are processed with the client's visit-func
    {visit} (if not NULL) and appended to the visit-list {vP} (if not
    NULL), according to the general rules. */

/* GENERAL ENUMERATION */

bool_t enum_items
  ( ref_vec_t root, 
    uint32_t ns, 
    enum_step_t *step[], 
    enum_visit_t *visit[], 
    ref_vec_t *vP
  );
  /* Enumerates all items that can be reached from the `root' items
    {root.e[0..root.ne-1]} by all possible sequences of calls to
    the step-funcs {step[0..ns-1]}.  
    
    For each new item found in the enumeration, including the starting
    ones, the procedure calls one of the visit-funcs {visit[k](p)},
    where {k} is some integer in {0..ns}. If {vP} is not NULL, the
    procedure will store into {*vP} all visited elements. See the
    comments above for more details on the parameters {step}, {visit},
    and {vP}.
    
    If any entry {visit[k]} is NULL, it is treated as a procedure that
    does nothing and returns FALSE.
    
    The enumeration stops when a call to any visit-func {visit[k]}
    returns TRUE, or when all reachable items have been visited. The
    {enum_orbits} procedure itself returns TRUE in the first case, and
    FALSE otherwise.
    
    Note that the {step[]} vector should have {ns} elements, while the
    {visit[]} vector should have {ns+1}. Whenever {visit[k](p)} is
    called with {k < ns}, then {step[k]} was the last step along the
    path that led the procedure to {p} from the {root} items. If
    {k==ns}, it means that the path is empty --- that is, {p} is a
    {root} item that is not reachable from the earlier ones through
    {step} chains. The first call to {visit}, in particular, always
    has {k==ns}.
    
    The enumeration is `stratified depth-first', with {step[i]}
    having higher enumeration priority than {step[j]} when {0 <= i < j
    <= ns-1}. More precisely: whenever {visit[k](p)} is called, the
    procedure will have visited every item that is reachable from a
    previously visited item by chains of steps of level less than {k}.
    After each call {visit[k](p)}, all items that are still unvisited
    and are are reachable from {p} through chains of steps with levels
    less than {k} will be visited next, before the next {visit[j]}
    call with {j >= k}.
    
    In particular, if {step[k]} is one-to-one, for all {k} in
    {0..ns-1}, then a call {visit[k](p)} means that {p} is the
    first item in a new orbit of the group generated by {{step[j] : j
    in 0..k-1 }}. moreover, that entire orbit will be visited next,
    before any subsequent {visit[j]} with {j>=k}. */

#endif
