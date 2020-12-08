#ifndef riso_H
#define riso_H

#define riso_DESC "Isomorphism-reduced labeled directed acyclic graphs"

/* Last edited on 2009-11-03 13:40:00 by stolfi */
/*©*/

#include <stdint.h>

#include <bool.h>
#include <vec.h>

/* 
  LETTERS, WORDS AND LANGUAGES: A /letter/ is a natural number. A
  word is a finite sequence letters. A /language/ is a finite set of
  words. */

typedef uint32_t riso_symbol_t;
  /* An arbitrary letter. */

/*
  ALPHABET: If {L} is a language, we denote by {\Sig L} its
  /alphabet/, the set of all letters that occur in the words of {L}.

  VOID and UNIT LANGUAGES: The /void language/ is the language {{}}
  with no words.  The /unit language/ is the language {{()}} that has
  only the empty word {()}.

  WEIGHT: The /weight/ of a language {L} is the number of its words
  plus the total number of letters in all those words. Note that the
  void language has weight 0, the unit language has weight 1, and
  all other languages have weight greater than 1.
  
  RECODING: A /recoding/ for a set {X} of letters is an injective
  function {f} from {X} to the natural numbers; and therefore a
  bijection from {X} to the set {f(X)}.
  
  A recoding extends naturally to a bijection from words over {X} to
  words over {f(X)}, and from languages over {X} to languages over
  {f(X)}. For any letter {a}, word {w}, and language {L},
  we write {a*f}, {w*f}, {L*f} for {f(a)}, {f(w)}, and {f(L)},
  respectively. */
  
typedef struct recoding_t 
  { riso_symbol_t in;
    riso_symbol_t ot;
    struct recoding_t *rest;
  } recoding_t;
  /* A recoding is represented by a linked list of {redoding_t}
    nodes, sorted by strictly decreasing order of their {in}-letters.
    The recoding maps {in} to {ot}, and letters less than {in} as specified 
    by {rest}. A {NULL} recoding stands for the unique recoding with
    empty domain.
    
    !!!{Make into a node of a larger structure} */

riso_symbol_t recoding_map_letter(recoding_t *rc, riso_symbol_t in);
  /* Maps the letter {in} according to the specified recoding. */
  
recoding_t *recoding_pair_set(recoding_t *rc, riso_symbol_t in, riso_symbol_t ot);  
   /* Returns a recoding {rc} that maps {in} to {ot}, and 
     other 

/*
  LANGUAGE ISOMORPHISM: Two languages {K,L} are said to be /isomorphic/ iff
  there is a recoding from {\Sig K} to {\Sig L} that extends to a
  bijection of {K} to {L}. Such a recoding is called and
  /isomorphism/ from {K} to {L}.

  Note that there may be two or more distinct isomorphisms between
  two languages, and that a language may be isomorphic to itself
  under an isomorphism that is not the identity (an /automorphism/).
  For example, the language {{00,11}} over the digits {0,1} is
  isomorphic to itself under the recoding {0-->1,1-->0}. */

/*
  PLEX: A /plex/ is a data structure that represents a 
  collection of languages over some alphabet.  */
 
typedef struct riso_t riso_t;
  /* A {riso_t} is a mutable plex. */
 
/*
  LEXICON: A /lexicon/ {P} is a node of a plex that represents a
  finite language, denoted by {<P>}. We write {\Sig P} for the
  alphabet {\Sig <P>} of that language.

  LEXICON EQUIVALENCE: Two lexicons {P,Q} are said to be /language
  equivalent/ iff they represent the same language; that is, iff
  {<P>==<Q>}.

  KINDS OF LEXICONS: We deal with three distict kinds of lexicons:
  /states/, /nodes/, and  /branches/. Some of them have 
  restrictions on the languages that they represent:

    STATES: The language of a state {S} may be the void language {{}}. 
    It may or may not have the empty word {()}.

    NODES: The language of a node {N} may be void, but may not include
    the empty word {()}. That is, every word has at least one
    letter.

    BRANCHES: The language of a branch {B} is not void, does not
    include the empty word, and all its words begin with the same
    letter --- the /anchor/ of {B}, denoted by {B.a}.

  VOID STATES AND NODES: A state or node {P} is /void/ iff {<P>}
  is the void language.  Note that a branch is never void.

  STRUCTURE OF A STATE: A state {P} is defined as a pair {(N,e)}
  where {N} is a node and {e} is a boolean. By definition, the
  language {<P>} is {<N>} if {e} is false, and {<N>\uni {()}} if {e}
  is true.
  
  STRUCTURE OF A NODE: A node {N} is defined as a set
  {{N.B[0],..N.B[m-1]}} of branches with distinct anchors. By
  definition, the language {<N>} is the (disjoint) union of the sets
  {<N.B[i]>} represented by its branches.

  Note that since branches are non-void, a node {N} is void iff it
  has no branches (that is, {m==0}).

  HEAD AND TAIL OF A NODE: The /head/ of a non-void node is its
  branch with highest anchor letter. Its /tail/ is the set of all
  other branches.

  STRUCTURE OF A BRANCH: A branch {B} is defined as a triple
  {(a,P,s)} where {a} is the anchor letter, {P} is a non-void state,
  and {s} is a recoding for {P}. The language {<B>} consist of all
  words in {<P>}, mapped by the recoding {s}, each prefixed with
  the letter {a}. */
   
typedef uint32_t riso_node_t;
  /* A node of a plex. */
  
typedef uint32_t riso_state_t;
  /* A state of a plex. */
     
typedef struct riso_branch_t 
  { riso_symbol_t in;
    riso_state_t st;
    riso_recoding_t *rc;
  } riso_branch_t;
  /* A branch of a plex. */

bool_t riso_state_is_unitary(riso_t *D, riso_state_t u);
riso_node_t riso_state_node(riso_t *D, riso_state_t u);
riso_state_t riso_state_find(riso_t *D, bool_t ac, riso_node_t nd);
     
bool_t riso_branch_anchor(riso_t *D, riso_branch_t br);
riso_state_t riso_branch_state(riso_t *D, riso_branch_t br);
riso_recoding_t riso_branch_recoding(riso_t *D, riso_branch_t br);
riso_branch_t riso_branch_find(riso_t *D, riso_symbol_t in, riso_state_t st, riso_recoding_t *rc);  

typedef struct riso_node_data_t 
  { riso_node_t tail;
    riso_branch_t head;
  } riso_node_data_t;
  /* The four attributes of a proper node. */

/*
  STRUCTURAL ORDERING: There is a weak total order {\leq}, the
  /structural ordering/, among lexicons of the same class. This
  ordering implies a /structural congruence/ relation, namely {P} is
  congruent to {Q} iff {P\leq Q} and {Q\leq P}. We describe this
  ordering in terms of a /structural comparison/ function {cmp(,)},
  with result in {{-1,0,+1}}, which is defined recursively as
  follows.

  * If {a} and {b} are letters, {cmp(a,b)} is the sign of {a-b}.

  * If {d} and {e} are booleans, {cmp(d,e)} is zero if {d==e},
    and {-1} if {d} is false and {e} is true.

  * If {s} and {t} are recodings, {cmp(s,t)} is zero if {s==t},
    and {-1} if {s} is less than {t} in lexicographic order.

  * If {P} and {Q} are both void states, then {cmp(P,Q) == 0}.

  * If {P} is a void state and {Q} is a non-void state, then {cmp(P,Q) == -1}.

  * If {P} and {Q} are both non-void states, then {cmp(P,Q)} is
    the first nonzero outcome in the sequence {cmp(P.N,Q.N)}
    and {cmp(P.e,Q.e)}; or zero if both outcomes are zero.

  * If {A} and {B} are branches, then {cmp(A,B)} the first nonzero
    outcome in the sequence {cmp(A.a,B.a)}, {cmp(A.P,B.P)}, and
    {cmp(A.s,B.s)}; or zero if all three outcomes are zero.

  * If {M} and {N} are nodes, then {cmp(M,N)} the first nonzero
    outcome in the sequence {cmp(M.head,N.head)} and {cmp(M.tail,N.tail)};
    or zero if both outcomes are zero.

  CONGRUENCE UNIQUENESS: For space efficiency, it is desirable to
  ensure that there are no distinct but congruent lexicons. To this
  end, we choose a unique representation for the void and unit
  states. Then, whenever a lexicon is to be constructed from its
  components, we use a hash table to check whether there is another
  lexicon with the same fields; if so we return that lexicon.

  LEXICON ISOMORPHISM: Two lexicons {P,Q} are said to be
  /isomorphic/ iff their languages {<P>}, {<Q>} are isomorphic.

  CONGRUENCE AND EQUIVALENCE:  Note that structural
  congruence between two lexicons {P,Q} (namely, {cmp(P,Q)==0})
  implies language equivalence ({<P>==<Q>}). However, the converse
  is not true: the same finite language may be represented by many
  lexicons that are not structurally congruent.

  In particular, if {B=(a,P,s)} is any branch, and {t} is any
  non-trivial isomorphism of some state {Q} to {P}, then {B} will
  the equivalent --- but not congruent! --- to the branch
  {(a,Q,t*s)}, where {t*s} denotes the composition of the recodings
  {t} (applied first) with {s} (applied last).

  Indeed, if {<P>} has any non-trivial automorphism {t}, then 
  {B=(a,P,s)} will be wquivalent, but not congruent, to the 
  branch {(a,P,t*s)}.

  ISOCANONIC LEXICONS: We can use the structural comparison to choose
  a unique representative for each isomorphism class of lexicons.
  Namely a lexicon {P} is /isocanonical/ if {P\leq Q} for 
  any other lexicon {Q} isomorphic to {P}.

  Two isomorphic lexicons which are both isocanonical must be 
  congruent to each other. Note that the alphabet
  of an isocanonical lexicon must be a set of consecutive
  integers, starting from 0.

  REDUCED LEXICONS: We now define the subclass of /reduced/ lexicons
  for which equivalence does imply congruence.

    * A state {P=(N,e)} is /reduced/ iff its node {N} is reduced.

    * A node is /reduced/ if its branches are all reduced. 
      Note that a void state is always reduced.

    * A branch {B=(a,P,s)} is /reduced/ iff the state {P} is reduced
      and isocanonical, and {s \leq h*s} for any automorphism {h} of {P}.

  EQUIVALENCE OF REDUCED LEXICONS: Theorem: if two reduced states
  (or nodes, or branches) are equivalent, they are congruent.

  Proof: by induction on the weight of {<P>\uni<Q>}. The only
  non-trivial case is for branches. If two branches {(a,P,s)} and
  {(b,Q,t)} are equivalent, then we must obviously have {a==b}, and
  {<P>*s == <Q>*t}. Therefore {<P> = <Q>*t*s^{-1}}, that is, {P} is
  isomorphic to {Q}. Since both branches are reduced, we must have
  {cmp(P,Q)\leq 0} and {cmp(Q,P)\leq 0}, that is, {cmp(Q,P) == 0}.

  ISOMORPHISMS OF REDUCED LEXICONS: To find the reduced lexicon for
  a given language, we need to be able to find the isocanonical
  state {Q} of a given reduced state {P}. The problem trivially
  reduces to finding the reduced isocanonical node {M} of a given
  node {N}.

  A void node is always isocanonical, so assume that {M}
  has {m>0} branches.  It is easy to see that the branch {M.head}
  must have {m-1} as the anchor letter. 

  Let {f} be the (unknown) isomorphism from {N} to {M}. Then each
  branch {A[i]=(a[i],P[i],s[i])} of {N} corresponds to a unique
  branch {B[j]=(b[j],Q[j],t[j])} of {M}, where {b[j] = a[i]*f} and
  the state {P[i]*s[i]*f} is equivalent to the state {Q[j]*t[j]}.
  Therefore {s[i]*f*t[j]^{-1}} is an isomorphism from {P[i]} to
  {Q[j]}. Since {A} and {B} are reduced, the states {P[i]} and
  {Q[j]} are both isocanonical; therefore {P[i]==Q[j]}.

  Therefore the states {Q[0],...,Q[m-1]} in the branches of node {M}
  are congruent (and not merely isomorphic) to the states
  {P[0],...P[m-1]} in the branches of node {N}. Because of the way
  nodes are compared by {cmp}, for {M} to be isocanonical the state
  {M.head.P} must be the minimum (in the sense of {cmp})
  among the states {P[0],...P[m-1]}.

  If that state is unique, then we just need to build the state
  {letters {b[0],...b[m-1]} must be the integers {0,..m-1}, in that
  order. 

  It remains to find the recodings {t[0],..t[m-1]}.  These must be the 
  lexicographically smallest recodings such that {<M>} is isomorphic to 
  {<N>} 
  
  A plex can be used to define a regular deterministic
  length-preserving mapping from strings to strings; or, in
  particular, to cassify a finite set of strings into finitely
  many classes; or simply to represent a finite set of strings. */

uint32_t riso_node_count (riso_t *D);
  /* Number of nodes currently in the plex {D}.
    Cost: {O(1)} time, 0 space. */

riso_node_t riso_node_max(riso_t *D);
  /* The current nodes of {D} are numbered from 0 to {riso_node_max(D)}.
    Cost: {O(1)} time, 0 space. */
 
/*
  NULL NODE AND PROPER NODES: Every plex has a /null node/, which is
  node number 0. Nodes that are not null are /proper/. */ 
 
#define riso_node_NULL (0)
 
/* 
  NODE LABELS: every proper node {s} has two labels: an /input mark/
  or /i-mark/, and an /output mark/ or /o-mark/. Each mark is an
  unsigned integer between 0 and some maximum value, that depends on
  the plex and may be different for i-marks and o-marks. */ 
  
riso_symbol_t riso_in(riso_t *D, riso_node_t s);
riso_symbol_t riso_ot(riso_t *D, riso_node_t s);
  /* The input and output marks of a node {s}.
    Cost: {O(1)} time, 0 space. */

riso_symbol_t riso_in_max(riso_t *D);
riso_symbol_t riso_ot_max(riso_t *D);
  /* The maximum valid input and output marks of a node of {D}.
    Cost: {O(1)} time, 0 space. */

/* 
  RELATED NODES: Every proper node {s} also has two associated nodes,
  the /links/ of {S}: the /fail link/ or /f-link/ and the /pass link/
  or /p-link/. */ 
 
riso_node_t riso_p_link (riso_t *D, riso_node_t s);
riso_node_t riso_tail (riso_t *D, riso_node_t s);
  /* The pass and fail links of a node {s}.
    Cost: {O(1)} time, 0 space. */
  
void riso_node_data_get(riso_t *D, riso_node_t s, riso_node_data_t *dt);  
  /* Stores into {*dt} the attributes of the proper node {s} of {D}. */

riso_node_t riso_node_find(riso_t *D, riso_node_data_t *dt);
  /* Returns a node in {D} with the node data {*dt}. If 
    there is no such node, return {riso_node_SKIP}.  Rebuilds 
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
  we always have {riso_in(D,s) > riso_in(D,t)}.
  
  NODE UNIQUENESS: The implementation ensures that nodes are uniquely
  represented, in the sense that no two distinct nodes have the same
  marks and links. That is, the implementation ensures the following
  invariant:
  
  { 
    s != t  ==> 
      \/ (s == 0) != (t == 0)
      \/ riso_in(D,s) != riso_in(D,t) 
      \/ riso_ot(D,s) != riso_ot(D,t) 
      \/ riso_p_link(D,s) != riso_tail(D,t) 
      \/ riso_tail(D,s) != riso_tail(D,t)
  } */
 
/* 
  SUBNODES: by definition, a /subnode/ of a node {s} of the plex is a
  node {t} that can be obtained from {s} by following zero or more f-links.
  Note that a node has only finitely many subnodes, and that 
  the null node is a subnode of every node. */

riso_node_t riso_subnode_find(riso_t *D, riso_node_t s, riso_symbol_t i);
  /* Follows f-links from {s} until it finds a proper subnode that 
    has input mark equal to {i}.  If there is no such subnode
    (in particular, if {s} is NULL), returns {riso_node_NULL}.
    Cost: {O(m)} time, 0 space, where {m} is the number of
    subnodes of {s}. */

uint32_t riso_subnode_count(riso_t *D, riso_node_t s);
  /* Number of arcs out of the automaton state represented by
    node {s} of {D}. It is zero iff {s} is the null node.
    Cost: {O(m)} time, 0 space, where {m} is the number of
    subnodes of {s}. */
 
riso_node_t riso_subnode_first(riso_t *D, riso_node_t s);
  /* The unique subnode of {s} that has only one arc out of it (which
    is the first arc out of the state represented by {s}). Fails if
    {s} is the null node. Cost: {O(m} time, 0 space, where {m} is the
    number of subnodes of {s}. */

/* 
  !!! Rename variables {t<->s} !!! 
  PATHS: A /path/ in a plex {D} starts with an output symbol {o[0]} and
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
  { riso_disp_FINE,
    riso_disp_SKIP,
    riso_disp_STOP
  } riso_disp_t;
  /* Return codes for {riso_enum_paths} and its argument actions. */

typedef riso_disp_t riso_node_action_t(uint32_t len, riso_symbol_t o, riso_node_t t);
  /* Type of a procedure that is called by {riso_enum_paths} when the
    current path has length {len} and ends with an o-symbol {o} at a
    node {t}. */

typedef riso_disp_t riso_p_step_action_t
  ( uint32_t len,
    riso_node_t s,
    riso_symbol_t i, 
    riso_symbol_t o,
    riso_node_t t
  );
  /* Type of a procedure that is called by {riso_enum_paths} when the
    current path has length {len} and is about to be extended 
    with the p-step {(s,i,o,t)}. */
  
typedef riso_disp_t riso_string_action_t(uint32_t nstr, riso_symbol_t str[]); 
  /* Type of a procedure that operates on a string {str[0..nstr-1]}
    of symbols. */

vec_typedef(riso_symbol_vec_t,riso_symbol_vec,riso_symbol_t);

riso_disp_t riso_enum_paths
  ( riso_t *D, 
    riso_symbol_t o_root,
    riso_node_t t_root,
    riso_node_action_t *enter,
    riso_p_step_action_t *push, 
    riso_p_step_action_t *pop,
    riso_node_action_t *exit
  );
  /*
    Enumerates a set of paths in the plex, starting from node {t_root},
    in double depth-first order.  The set of paths is defined by the
    client-provided procedures {enter}, {push}, {pop}, and {exit}.
    Only the p-link steps of the path are reported, the f-link
    steps being left implied.

    The enumeration algorithm can be described in terms of a
    conceptual {current path} that grows or shrinks one p-step at a time.
    The path initially has length 0 and contains just the symbol {o[0]
    = o_root} and the node {t[0] = t_root}. The procedure
    {riso_enum_paths} will call:

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

    The action procedures should normally return {riso_disp_FINE}.
    If any of them returns {riso_disp_STOP}, the enumeration is 
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
    
    In particular, {riso_enum_paths} will call {enter(0,o_root,t_root)} at the very
    beginning, and (if not aborted) will call {exit(0,o_root,t_root)} at the very end.

    Note that the same node {t} may be {enter}ed and {exit}ed many, many times.

    The client procedures can prune branches of the path tree by returning
    {riso_disp_SKIP}.  Specifically,

      if {enter(len,o,t)} returns {riso_disp_SKIP}, then
        {riso_enum_paths} will omit the entire branch of the path tree
        rooted at current path, and call {exit(len,o,t)} right away,
        as if {t} had no outgoing arcs;

      if {push(len,s,i,o,t)} returns {riso_disp_SKIP}, then
        {riso_enum_paths} will call {pop(len,s,i,o,t)} right away;

      if {pop(len,s,i,o,t)} returns {riso_disp_SKIP}, then
        {riso_enum_paths} ignores any remaining arcs out of the 
        last entered node {t'}, and calls {exit(len,o',t')} right away.
      
      if {exit(len,o,t)} returns {riso_disp_SKIP}, the effect
        is the same as if it had returned {riso_disp_FINE}.

    If any action procedure is NULL, {riso_enum_paths} will provide a
    a trivial default action that does nothing and returns {riso_disp_FINE}. */

/* 
  CONSTRUCTION: The following procedures may be used build an automaton
  with all invariants satisfied. */

#define riso_ni_MAX (32)
#define riso_no_MAX (32)
#define riso_nn_MAX (30)
#define riso_ntot_MAX (64)

riso_t *riso_new(uint32_t nn, uint32_t ni, uint32_t no, riso_node_t max_alloc_node);
  /* Creates a new {riso_t}, with node numbers of {nn} bits, input symbols
    of {ni} bits, and output symbols of {no} bits. These parameters
    must not exceed {riso_nn_MAX}, {riso_ni_MAX}, and {riso_no_MAX},
    respectively; and their sum must not exceed {riso_ntot_MAX}.

    The plex will initially have no proper nodes, only
    {riso_node_NULL}. The parameter {max_alloc_node} is the maximum
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

riso_node_t riso_node_from_fields(riso_t *D, riso_node_t f, riso_symbol_t i, riso_symbol_t o, riso_node_t p);
  /* Returns the unique node {s} with f-link {f}, i-mark {i}, o-mark
    {o}, and p-link {p}. Cost: {O(1)} time and space, apart from
    eventual expansion and rehasing costs. */

riso_node_t riso_node_from_data(riso_t *D, riso_node_data_t *dt);
  /* Returns the unique node {s} with the node data {*dt} Equivalent
    to {riso_node_from_fields(D,dt.f,dt.i,dt.o,dt.p)}. Cost: {O(1)}
    time and space, apart from eventual expansion and rehasing
    costs. */

riso_node_t *riso_copy(riso_t *OLD, riso_t *NEW, riso_node_t old_s, riso_node_t map[]);
  /*
    Copies into the {NEW} plex all nodes of the {OLD} plex that are
    reachable from {old_s}. 

    If the {map} argument is not NULL, it is assumed to be a table
    that tells which nodes of {OLD} have been copied to {NEW}, and
    which are their corresponding nodes in {NEW}. Namely, {map} must
    have at least {riso_node_max(OLD)} elements. For any proper node
    {s} of {OLD}, if {map[s-1]} is positive, the procedure assumes
    that node {s} of {OLD} is equivalent to node {map[s-1]} of {NEW}.
    In that case the p-link and the f-link of {s} in {OLD} must have
    been copied also, and their {map} entries must reflect that.
    
    If {map} is NULL, {riso_copy} will allocate a vector of
    the appropriate size. 
    
    In any case, {riso_copy} will set {map[s-1]} appropriately for every
    node {s} of {OLD} that it copies, and will return the address of the map
    (whether user-provided or allocated internally).

    The procedure will maintain the uniqueneness invariant: when
    copying a node, it will create a new node in {NEW} only if there
    is no equivalent node there. This will be true even if the
    existence of that node was not recorded in the client-given
    {map}. */

void riso_discard(riso_t *D, riso_node_t s);
  /* Discards node {s} and all higher-numbered nodes. Invalidates the
    hash table. */ 


/* STORAGE MANAGEMENT */

uint32_t riso_node_max_alloc(riso_t *D);
  /* Max node for which storage has been allocated. Clients can create
    at least {riso_node_max_alloc - riso_node_max} new nodes
    without triggering an automatic expansion. */
 
void riso_expand(riso_t *D, riso_node_t max_alloc_node);
  /* Expands the internal storage to accomodate nodes
    from 0 to {riso_node_max_alloc}. Preserves all current nodes
    and their numbers. Bombs out if there is no more space,
    or if {riso_node_max_alloc} exceeds the maxmum possible 
    node number. */

uint32_t riso_reachable_node_count(riso_t *D, uint32_t nroots, riso_node_t root[]);
  /* Count of all distinct proper nodes that are reachable from nodes
    {root[0..nroots-1]} by any sequence of f-links or p-links. The
    count includes the roots themselves, excluding repeated or null
    nodes; and is zero if and only if all roots are null. Cost: {O(m)}
    time and space, where {m} is the maximum root node. */
 
void riso_crunch(riso_t *D, uint32_t nroots, riso_node_t root[]);
  /* Discards all nodes not reachable from the nodes {root[0..nroots-1]},
    squeezes the reachable ones together (preserving their order), and
    updates the {root} vector to reflect the new node numbering.

    WARNING: By definition, {riso_crunch} generally changes the
    numbers of all nodes, and deletes some of them. The client must make
    sure that all important node variables are part of the {root}
    vector. */
    
void riso_free(riso_t *D);
  /* frees all storage used by {D}, including the record {*D}
    itself and all internal tables. */
    
/* UNSAFE PROCEDURES */

void riso_node_data_set(riso_t *D, riso_node_t s, riso_node_data_t *dt);  
  /* Stores the attributes {*dt} for the node {s}, which must be
    proper. Checks field value ranges, topological ordering and
    lexicographic ordering; but does not check uniqueness, an
    invalidates the hash table. */
  
#endif
