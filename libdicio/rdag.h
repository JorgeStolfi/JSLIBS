#ifndef rdag_H
#define rdag_H

#define rdag_DESC "Reduced labeled directed acyclic graphs"

/* Last edited on 2009-10-30 22:27:43 by stolfi */
/*©*/

#include <stdint.h>

#include <bool.h>
#include <vec.h>

/* 
  DAG: A /reduced labeled directed acyclic graph/, here called simply
  /dag/, is an entity equivalent to a reduced deterministic acyclic
  finite-state semiautomaton.
  
  A dag can be used to define a regular deterministic
  length-preserving mapping from strings to strings; or, in
  particular, to cassify a finite set of strings into finitely
  many classes; or simply to represent a finite set of strings. */
  
typedef struct rdag_t rdag_t;
  /* An {rdag_t} is a mutable reduced labeled dag. */

/* NODES: A dag consists of a finite set of /nodes/. In a dag with {N}
  nodes, they are the unsigned integers in the range {0} to {N-1}. */
  
typedef uint32_t rdag_node_t;
  /* A node of a dag. */

uint32_t rdag_node_count (rdag_t *D);
  /* Number of nodes currently in the dag {D}.
    Cost: {O(1)} time, 0 space. */

rdag_node_t rdag_node_max(rdag_t *D);
  /* The current nodes of {D} are numbered from 0 to {rdag_node_max(D)}.
    Cost: {O(1)} time, 0 space. */
 
/*
  NULL NODE AND PROPER NODES: Every dag has a /null node/, which is node number 0.
  Nodes that are not null are /proper/. */ 
 
#define rdag_node_NULL (0)
 
/* 
  NODE LABELS: every proper node {s} has two labels: an /input mark/
  or /i-mark/, and an /output mark/ or /o-mark/. Each mark is an
  unsigned integer between 0 and some maximum value, that depends on
  the dag and may be different for i-marks and o-marks. */ 

typedef uint32_t rdag_symbol_t;
  /* An arbitrary symbol. */
  
rdag_symbol_t rdag_i_mark(rdag_t *D, rdag_node_t s);
rdag_symbol_t rdag_o_mark(rdag_t *D, rdag_node_t s);
  /* The input and output marks of a node {s}.
    Cost: {O(1)} time, 0 space. */

rdag_symbol_t rdag_i_mark_max(rdag_t *D);
rdag_symbol_t rdag_o_mark_max(rdag_t *D);
  /* The maximum valid input and output marks of a node of {D}.
    Cost: {O(1)} time, 0 space. */

/* 
  RELATED NODES: Every proper node {s} also has two associated nodes,
  the /links/ of {S}: the /fail link/ or /f-link/ and the /pass link/
  or /p-link/. */ 
 
rdag_node_t rdag_p_link (rdag_t *D, rdag_node_t s);
rdag_node_t rdag_f_link (rdag_t *D, rdag_node_t s);
  /* The pass and fail links of a node {s}.
    Cost: {O(1)} time, 0 space. */
   
typedef struct rdag_node_data_t 
  { rdag_node_t f_link;
    rdag_symbol_t i_mark;
    rdag_symbol_t o_mark;
    rdag_node_t p_link;
  } rdag_node_data_t;
  /* The four attributes of a proper node. */
  
void rdag_node_data_get(rdag_t *D, rdag_node_t s, rdag_node_data_t *dt);  
  /* Stores into {*dt} the attributes of the proper node {s} of {D}. */

rdag_node_t rdag_node_find(rdag_t *D, rdag_node_data_t *dt);
  /* Returns a node in {D} with the node data {*dt}. If 
    there is no such node, return {rdag_node_SKIP}.  Rebuilds 
    the hash tables if they are not valid.
    
    Cost: if the hash tables are valid, {O(1)} time and space. If the
    hash tables are not valid, {O(N)} time and space, where {N}
    is the number of nodes in {D}. */
 
/*
  TOPOLOGICAL ORDERING: For every proper node {s}, the pass
  node {p} and the fail link {f} are strictly less than {s}. It
  follows that from every node {s} one will eventually reach the null
  node after following a finite number of f- or p-links, no matter
  how they are interleaved.
    
  LEXICOGRAPHIC ORDERING: For every proper node {s} with a proper f-link {t}
  we always have {rdag_i_mark(D,s) > rdag_i_mark(D,t)}.
  
  NODE UNIQUENESS: The implementation ensures that nodes are uniquely
  represented, in the sense that no two distinct nodes have the same
  marks and links. That is, the implementation ensures the following
  invariant:
  
  { 
    s != t  ==> 
      \/ (s == 0) != (t == 0)
      \/ rdag_i_mark(D,s) != rdag_i_mark(D,t) 
      \/ rdag_o_mark(D,s) != rdag_o_mark(D,t) 
      \/ rdag_p_link(D,s) != rdag_f_link(D,t) 
      \/ rdag_f_link(D,s) != rdag_f_link(D,t)
  } */
 
/* 
  SUBNODES: by definition, a /subnode/ of a node {s} of the dag is a
  node {t} that can be obtained from {s} by following zero or more f-links.
  Note that a node has only finitely many subnodes, and that 
  the null node is a subnode of every node. */

rdag_node_t rdag_subnode_find(rdag_t *D, rdag_node_t s, rdag_symbol_t i);
  /* Follows f-links from {s} until it finds a proper subnode that 
    has input mark equal to {i}.  If there is no such subnode
    (in particular, if {s} is NULL), returns {rdag_node_NULL}.
    Cost: {O(m)} time, 0 space, where {m} is the number of
    subnodes of {s}. */

uint32_t rdag_subnode_count(rdag_t *D, rdag_node_t s);
  /* Number of arcs out of the automaton state represented by
    node {s} of {D}. It is zero iff {s} is the null node.
    Cost: {O(m)} time, 0 space, where {m} is the number of
    subnodes of {s}. */
 
rdag_node_t rdag_subnode_first(rdag_t *D, rdag_node_t s);
  /* The unique subnode of {s} that has only one arc out of it (which
    is the first arc out of the state represented by {s}). Fails if
    {s} is the null node. Cost: {O(m} time, 0 space, where {m} is the
    number of subnodes of {s}. */

/* 
  !!! Rename variables {t<->s} !!! 
  PATHS: A /path/ in a dag {D} starts with an output symbol {o[0]} and
  a node {t[0]}, and continues with zero or more /p-steps/. Each
  p-step is a quadruple {(s[k],i[k],o[k],t[k])} of {D}, for {k} in
  some range {1..n}; such that each {s[k]} is a subnode of {t[k-1]}
  with i-mark {i[k]}, o-mark {o[k]}, and p-link {t[k]}. The integer
  {n} is the /length/ of the path.
  
  STANDARD PATH ORDERS: The ordering of the subnodes of each node
  induces an ordering of the set of all paths that begin with a given
  node {t}; namely, lexicographic in the sequence of input symbols
  {i[1],i[2],...,i[n]}.
  
  There are actually three lexicographic path orderings, differing on
  whether the paths that /end/ at {o,t} appear either before all the
  paths that /go through/ {o,t} (/pre-order/), or after them
  (/end-order/), or both before and after them (/double order/). */

typedef enum 
  { rdag_disp_FINE,
    rdag_disp_SKIP,
    rdag_disp_STOP
  } rdag_disp_t;
  /* Return codes for {rdag_enum_paths} and its argument actions. */

typedef rdag_disp_t rdag_node_action_t(uint32_t len, rdag_symbol_t o, rdag_node_t t);
  /* Type of a procedure that is called by {rdag_enum_paths} when the
    current path has length {len} and ends with an o-symbol {o} at a
    node {t}. */

typedef rdag_disp_t rdag_p_step_action_t
  ( uint32_t len,
    rdag_node_t s,
    rdag_symbol_t i, 
    rdag_symbol_t o,
    rdag_node_t t
  );
  /* Type of a procedure that is called by {rdag_enum_paths} when the
    current path has length {len} and is about to be extended 
    with the p-step {(s,i,o,t)}. */
  
typedef rdag_disp_t rdag_string_action_t(uint32_t nstr, rdag_symbol_t str[]); 
  /* Type of a procedure that operates on a string {str[0..nstr-1]}
    of symbols. */

vec_typedef(rdag_symbol_vec_t,rdag_symbol_vec,rdag_symbol_t);

rdag_disp_t rdag_enum_paths
  ( rdag_t *D, 
    rdag_symbol_t o_root,
    rdag_node_t t_root,
    rdag_node_action_t *enter,
    rdag_p_step_action_t *push, 
    rdag_p_step_action_t *pop,
    rdag_node_action_t *exit
  );
  /*
    Enumerates a set of paths in the dag, starting from node {t_root},
    in double depth-first order.  The set of paths is defined by the
    client-provided procedures {enter}, {push}, {pop}, and {exit}.
    Only the p-link steps of the path are reported, the f-link
    steps being left implied.

    The enumeration algorithm can be described in terms of a
    conceptual {current path} that grows or shrinks one p-step at a time.
    The path initially has length 0 and contains just the symbol {o[0]
    = o_root} and the node {t[0] = t_root}. The procedure
    {rdag_enum_paths} will call:

        {enter(len,o,t)} whenever the current path, which is being
            enumerated for the first time, has length {len} and ends
            with output mark {o} and node {t}.

        {push(len,s,i,o,t)} whenever the current path has {len} arcs,
            and is about to be extended with a p-step {(s,i,o,t)}.

        {pop(len,s,i,o,t)} whenever the current path has length {len} arcs
            and has just lost its final p-step {(s,i,o,t)}.

        {exit(len,o,t)} whenever the current path, which is being
            enumerated for the last time, has length {len} and
            ends with output symbol {o} and node {t}.

    The action procedures should normally return {rdag_disp_FINE}.
    If any of them returns {rdag_disp_STOP}, the enumeration is 
    terminated and the procedure returns immediately, withou calling
    any other action.

    Unless the enumeration is stopped, every call to {enter} will be followed
    eventually by a matching call to {exit}, and every call to {push}
    will be followed eventually by a matching call to {pop}.  The typical
    action pattern for a generic node {s} with {n} outgoing arcs is

    |     enter(len,o,t)
    |       push(len,s[1],i[1],o[1],t[1])
    |         ...
    |       pop(len,s[1],i[1],o[1],t[1])
    |       push(len,s[2],i[2],o[2],t[2])
    |         ...
    |       pop(len,s[2],i[2],o[2],t[2])
    |       ...
    |       push(len,s[k],i[k],o[k],t[k])
    |         ...
    |       pop(len,s[k],i[k],o[k],t[k])
    |     exit(len,o,t)
    
    where {(i[1],o[1],t[1]), (i[2],o[2],t[2]), ..., (i[k],o[k],t[k])}
    are the arcs out of state {t}, and {s[1],s[2],...,s[k]} are the 
    corresponding subnodes of {t}. 
    
    In particular, {rdag_enum_paths} will call {enter(0,o_root,t_root)} at the very
    beginning, and (if not aborted) will call {exit(0,o_root,t_root)} at the very end.

    Note that the same node {t} may be {enter}ed and {exit}ed many, many times.

    The client procedures can prune branches of the path tree by returning
    {rdag_disp_SKIP}.  Specifically,

      if {enter(len,o,t)} returns {rdag_disp_SKIP}, then
        {rdag_enum_paths} will omit the entire branch of the path tree
        rooted at current path, and call {exit(len,o,t)} right away,
        as if {t} had no outgoing arcs;

      if {push(len,s,i,o,t)} returns {rdag_disp_SKIP}, then
        {rdag_enum_paths} will call {pop(len,s,i,o,t)} right away;

      if {pop(len,s,i,o,t)} returns {rdag_disp_SKIP}, then
        {rdag_enum_paths} ignores any remaining arcs out of the 
        last entered node {t'}, and calls {exit(len,o',t')} right away.
      
      if {exit(len,o,t)} returns {rdag_disp_SKIP}, the effect
        is the same as if it had returned {rdag_disp_FINE}.

    If any action procedure is NULL, {rdag_enum_paths} will provide a
    a trivial default action that does nothing and returns {rdag_disp_FINE}. */

/* 
  CONSTRUCTION: The following procedures may be used build an automaton
  with all invariants satisfied. */

#define rdag_ni_MAX (32)
#define rdag_no_MAX (32)
#define rdag_nn_MAX (30)
#define rdag_ntot_MAX (64)

rdag_t *rdag_new(uint32_t nn, uint32_t ni, uint32_t no, rdag_node_t max_alloc_node);
  /* Creates a new {rdag_t}, with node numbers of {nn} bits, input symbols
    of {ni} bits, and output symbols of {no} bits. These parameters
    must not exceed {rdag_nn_MAX}, {rdag_ni_MAX}, and {rdag_no_MAX},
    respectively; and their sum must not exceed {rdag_ntot_MAX}.

    The dag will initially have no proper nodes, only
    {rdag_node_NULL}. The parameter {max_alloc_node} is the maximum
    node number that the automaton may come to have before an
    automatic expansion. */

/* 
  AUGMENTATION: The following procedures add new nodes to an
  automaton. Existing nodes are not affected. Occasionally, they may
  cause some internal tables to be automatically expanded, and hash
  tables to be rebuilt. Since the expansion increases the capacity
  roughly in geometric progression, the amortized cost of expansion
  and rehashing is usually {O(1)} per operation, in the amortized
  sense. */

rdag_node_t rdag_node_from_fields(rdag_t *D, rdag_node_t f, rdag_symbol_t i, rdag_symbol_t o, rdag_node_t p);
  /* Returns the unique node {s} with f-link {f}, i-mark {i}, o-mark
    {o}, and p-link {p}. Cost: {O(1)} time and space, apart from
    eventual expansion and rehasing costs. */

rdag_node_t rdag_node_from_data(rdag_t *D, rdag_node_data_t *dt);
  /* Returns the unique node {s} with the node data {*dt} Equivalent
    to {rdag_node_from_fields(D,dt.f,dt.i,dt.o,dt.p)}. Cost: {O(1)}
    time and space, apart from eventual expansion and rehasing
    costs. */

rdag_node_t *rdag_copy(rdag_t *OLD, rdag_t *NEW, rdag_node_t old_s, rdag_node_t map[]);
  /*
    Copies into the {NEW} dag all nodes of the {OLD} dag that are
    reachable from {old_s}. 

    If the {map} argument is not NULL, it is assumed to be a table
    that tells which nodes of {OLD} have been copied to {NEW}, and
    which are their corresponding nodes in {NEW}. Namely, {map} must
    have at least {rdag_node_max(OLD)} elements. For any proper node
    {s} of {OLD}, if {map[s-1]} is positive, the procedure assumes
    that node {s} of {OLD} is equivalent to node {map[s-1]} of {NEW}.
    In that case the p-link and the f-link of {s} in {OLD} must have
    been copied also, and their {map} entries must reflect that.
    
    If {map} is NULL, {rdag_copy} will allocate a vector of
    the appropriate size. 
    
    In any case, {rdag_copy} will set {map[s-1]} appropriately for every
    node {s} of {OLD} that it copies, and will return the address of the map
    (whether user-provided or allocated internally).

    The procedure will maintain the uniqueneness invariant: when
    copying a node, it will create a new node in {NEW} only if there
    is no equivalent node there. This will be true even if the
    existence of that node was not recorded in the client-given
    {map}. */

void rdag_discard(rdag_t *D, rdag_node_t s);
  /* Discards node {s} and all higher-numbered nodes. Invalidates the
    hash table. */ 


/* STORAGE MANAGEMENT */

uint32_t rdag_node_max_alloc(rdag_t *D);
  /* Max node for which storage has been allocated. Clients can create
    at least {rdag_node_max_alloc - rdag_node_max} new nodes
    without triggering an automatic expansion. */
 
void rdag_expand(rdag_t *D, rdag_node_t max_alloc_node);
  /* Expands the internal storage to accomodate nodes
    from 0 to {rdag_node_max_alloc}. Preserves all current nodes
    and their numbers. Bombs out if there is no more space,
    or if {rdag_node_max_alloc} exceeds the maxmum possible 
    node number. */

uint32_t rdag_reachable_node_count(rdag_t *D, uint32_t nroots, rdag_node_t root[]);
  /* Count of all distinct proper nodes that are reachable from nodes
    {root[0..nroots-1]} by any sequence of f-links or p-links. The
    count includes the roots themselves, excluding repeated or null
    nodes; and is zero if and only if all roots are null. Cost: {O(m)}
    time and space, where {m} is the maximum root node. */
 
void rdag_crunch(rdag_t *D, uint32_t nroots, rdag_node_t root[]);
  /* Discards all nodes not reachable from the nodes {root[0..nroots-1]},
    squeezes the reachable ones together (preserving their order), and
    updates the {root} vector to reflect the new node numbering.

    WARNING: By definition, {rdag_crunch} generally changes the
    numbers of all nodes, and deletes some of them. The client must make
    sure that all important node variables are part of the {root}
    vector. */
    
void rdag_free(rdag_t *D);
  /* frees all storage used by {D}, including the record {*D}
    itself and all internal tables. */
    
/* UNSAFE PROCEDURES */

void rdag_node_data_set(rdag_t *D, rdag_node_t s, rdag_node_data_t *dt);  
  /* Stores the attributes {*dt} for the node {s}, which must be
    proper. Checks field value ranges, topological ordering and
    lexicographic ordering; but does not check uniqueness, an
    invalidates the hash table. */
  
#endif
