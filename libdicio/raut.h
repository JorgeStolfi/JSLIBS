#ifndef raut_H
#define raut_H

#define raut_DESC "Reduced deterministic acyclic automata"

/* Last edited on 2009-10-31 00:00:12 by stolfi */
/*©*/

#include <stdint.h>

#include <bool.h>
#include <rdag.h>

/*
  ATOMATA: For this module, an /automaton/ is a mutable data structure
  that represents a deterministic acyclic finite automaton, reduced
  but not necessarily connected. */
  
typedef struct raut_t raut_t;
  /* A representation of an automaton. */

/*
  INPUT SYMBOLS: The input alphabet of an automaton is a set of integers
  from 0 to some maximum value.  Each symbol is an {rdag_symbol_t}.
  
  STATES: An automaton has a finite set of states, represented by
  pairs {(o,s)} where {o} is 0 or 1 (the state's /accept bit/) and {s}
  is a natural number (the state's /node/). At any moment, the valid
  nodes are a range from 0 to some maximum node number. Each node is
  an {rdag_node_t}.
  
  ACCEPTING AND NON-ACCEPTING STATES: A state {u = (o,s)} of {A} is said to be
  /accepting/ if its o-mark {o} is 1; and /non-accepting/
  of its o-mark is zero.
  
  VOID STATE AND UNIT STATE: In particular, every automaton has the
  /void state/ {(0,NULL)} and the /unit state/ {(1,NULL)}. */
  
typedef struct raut_state_t { rdag_symbol_t ac; rdag_node_t nd; } raut_state_t;
  /* A state of an automaton with node {nd}.  The field {ac} is either
    0 or 1. !!! Should use the encoding {(o,s) --> 2*s+o} !!!. */

rdag_node_t raut_node_max(raut_t *A);
  /* The maximum node number currently valid in the automaton {A}.
    Cost: {O(1)} time, 0 space. */

bool_t raut_state_is_accepting(raut_t *A, raut_state_t u);
  /* TRUE iff {u} is an accepting state of the automaton {A}.
    Cost: {O(1)} time, 0 space. */

#define raut_state_VOID ((raut_state_t){ .ac = 0, .nd = rdag_node_NULL })
  /* The unique non-accepting state with null node. */

#define raut_state_UNIT ((raut_state_t){ .ac = 1, .nd = rdag_node_NULL })
  /* The unique accepting state with null node. */
  
/*
  ROOT STATE: Every automaton also has a distinguished /root state/.
  It may be any state, including the void state or the unit state. */

raut_state_t raut_root_get(raut_t *A);

/*
  TRANSITIONS: Each state {u} of {A} has a list of /outgoing arcs/ (transitions).
  Each transition is a pair {(i,u')} where {i} is an input symbol and 
  {u'} is some other state of {A} other than the void state. In 
  this list, the input symbols {i} are all distinct and increasing.  
  
  The list of outgoing arcs of a state {u} is empty if and only if
  {u} is the void state or the unit state. */

bool_t raut_state_has_arcs(raut_t *A, raut_state_t u);
  /* TRUE iff the state {s} has at least one proper outgoing arc.
    Cost: {O(1)} time, 0 space. */
 
raut_state_t raut_arc_follow(raut_t *A, raut_state_t u, rdag_symbol_t i);
  /* If there is an arc {(i,v)} out of state {u}, for some state {v},
    returns that state {v} (which is never the void state). Otherwise
    returns the void state. */

/* 
  INPUT STRINGS: An /input string/ of the automaton {A} is a finite
  sequence of zero or more input symbols.
  
  PATHS: A path in the automaton {A} is a sequence {P=(u[0],a[1],u[1],...,a[n],u[n])}
  where each {u[k]} is a state and each {a[k]} is an arc out of state
  {u[k-1]} to state {u[k]}. The integer {n} (which may be zero) is the 
  /length/ of the path.
  
  SPELLING: A path {P=(u[0],a[1],u[1],...,a[n],u[n])} is said to /spell/ the string
  {(i[1],...,i[n]), and vice-versa. Note that an input string spells 
  at most one path in {A}. */

/*  
  ACCEPTED AND REJECTED STRINGS: An input string {x = (x[1],x[2],...,x[n])} is
  /accepted/ by a state {u} of {A} if there is a path in {A} that starts
  at {u}, spells the string {x}, and ends at some accepting state. 
  
  The string {x} is /rejected/ by state {u} if it is not accepted;
  that is, if an attempt to spell {x} on the automaton, starting from
  {u}, either ends at a non-accepting state, or reaches a `dead end'
  along the way: a state {v} which has no outgoing arc with the
  appropriate input symbol. 
  
  We also say that a string is /accepted/ or /rejected/ by the automaton
  if it is respectively accepted or rejected by its root state. */
 
bool_t raut_state_accepts(raut_t *A, raut_state_t u, uint32_t nstr, rdag_symbol_t str[]);
  /* Returns TRUE if and only iff the automaton {A}, staring at state {u},
    accepts the input string {str[0..nstr-1]}. */
  
/*
  SUFFIXES OF A STATE: each state {u} is said to {accept} or
  {recognize} a language {Suff(u)}, the set of /suffixes/ of {u},
  which is the set of all input strings that spell directed paths from
  {u} that end at an accepting state. In particular, {Suff(u)}
  contains the empty string if and only if {u} is an accepting state.
  Indeed, {Suff((0,NULL))} is the /empty language/ {{}} (with no
  strings), and {Suff((1,NULL))} is the /unit language/ {{()}}
  (whose only string is the empty string).

  STATE UNIQUENESS: the implementation ensures that the automaton
  is /reduced/, that is, no two distinct states (proper or null)
  have the same set of suffix strings.  
  
  In particular, the void state is the only state that accepts the
  void language {{}}; and the unit state is the only state that
  accepts the unit language {{()}}. */

rdag_disp_t raut_enum_suffs(raut_t *A, raut_state_t u, rdag_string_action_t *proc);
  /* Enumerates all strings in the suffix set {Suff(u)} of state {u} in {A}.
    Calls {proc(n,x)} for each string {x[0..n-1]}. 
    
    The strings are enumerated in lexicographic order; each string {x}
    is visited before all longer strings of {Suff(u)} that begin with {x}.
    If {proc} returns {rdag_disp_SKIP}, skips all suffixes that begin with
    {x[0..n-1]}. If {proc} returns {rdag_disp_STOP}, stops the enumeration
    immediately and returns {rdag_disp_STOP}.  In all other cases returns
    {rdag_disp_FINE}. */

/* 
  CREATING AN AUTOMATON */
  
raut_t *raut_new(uint32_t nn, uint32_t ni, rdag_node_t max_alloc_node);
  /* Creates a new reduced automaton, with node numbers of {nn} bits and
    input symbols of {ni} bits. These parameters must not exceed
    {rdag_nn_MAX} and {rdag_ni_MAX}, respectively; and their sum must
    not exceed {rdag_ntot_MAX - 1}

    The automaton will initially have only the void and unit states. The
    root is initially set to the void state. The parameter
    {max_alloc_node} is the maximum node number that the automaton may
    come to have before an automatic expansion. */

/* 
  MODIFYING AN AUTOMATON: These procedures modify the automaton by
  adding more states, without disturbing the existing ones. Their 
  cost is {O(nstr)} time and space if no expansion is neded; {O(N)} time and space
  if expandion is needed, where {N} is the number of states in {A};
  {O(1)} time and space in amortized sense. */

raut_state_t raut_arc_set(raut_t *A, raut_state_t u, rdag_symbol_t i, raut_state_t v); 
  /* Returns a state {w} of {A} which has the same arcs as {u} except that it goes
    to{v} on input {i}. More precisely:
    
    * If {v} is a non-void state, the new state {w} will
    have {(i,v)} as one of its outgoing arcs. In particular, if {u}
    already has that outgoing arc, the precedure will return {u}.
    
    * If {v} is the void state, then {w} will not have any outgoing
    arc with the input symbol {i}. In particular, if {u} already has
    that outgoing arc, the precedure will return {u}.
    
    For any other input symbol {i'}, either the states {w}
    and {u} will have the same arc {(i',v')}, with same destination
    {v'}; or both states will lack such arc.
    
    If the new state {w} is not currently in {A}, the proedure will
    add it (and possibly some of its substates), expanding {A} as
    needed. Apart from expansion, the cost is {O(m)} time and space,
    where {m} is the number of arcs out of {u} that have input symbols
    greater than or equal to {i}. */

raut_state_t raut_string_add(raut_t *A, raut_state_t u, uint32_t nstr, rdag_symbol_t str[]); 
  /* Returns a state of {A} whose suffix set is {Suff(u)} plus the input string
    {str[0..nstr-1]}.  Creates the state and expands the internal tables if needed. 
    Existing states are not modified.  In particular, if the string is already in 
    {Suff(u)}, returns {u} itself. */

raut_state_t raut_string_remove(raut_t *A, raut_state_t u, uint32_t nstr, rdag_symbol_t str[]);
  /* Returns a state of {A} whose suffix set is {Suff(u)} minus the input string
    {str[0..nstr-1]}.  Creates the state and expands the internal tables if needed. 
    Existing states are not modified.  In particular, if the string is not in 
    {Suff(u)}, returns {u} itself. */

void raut_root_set(raut_t *A, raut_state_t u);
  /* Makes {u} the current root state of {A}.  Cost: {O(1)} time, 0 space. */

/*
  VIRTUAL LINKS: When a state {u} of {A} lacks an output arc with a given input mark
  {i}, it is convenient to assume that it has a /virtual arc/ {(i,(0,NULL))},
  with input mark {i}, output mark 0, and p-link null.
  Moreover, it is convenient to assume that the null node itself 
  has a virtual output arc of this form for every possible i-mark {i}.
  With these assumptions, every string spells a path in {D};
  and a string {x} is rejected by {A} if and only if a path 
  that spells {x} from {u_root} ends at a non-accepting state. */
  
/*
  PREFIXES OF A STATE: each proper state {u} has also a set of
  /prefixes/ {Pref(u)}, the set of all input strings spelled by paths
  that lead from the root state {u_root} to {u}. In particular, if
  {Pref(u_root)} is the unit language {{()}}; and if {u} is not
  reachable from the root state, then {Pref(u)} is the empty language.
  
  For technical and logical reasons, the set {Pref((0,NULL))} is undefined.

  DETERMINISM: according to the definitions above, the outgoung arcs
  of any given state are labeled with distinct symbols.  It follows
  that the automaton is deterministic: given any state {u} and any
  String {w}, there is exactly one state {t} (possibly null) that
  can be reached from {u} by a path whole labels spell {w}.
  Therefore, any two distinct states have disjoint prefix sets;
  and any proper reachable state can be uniquely identified by any
  of its prefix strings.
  
  SUBSTATES: A {substate} of a state {u=(o,s)} of {A} is any
  state {v=(o,t)} whose o-mark {o} is the same as that of {u},
  and whose node {t} is a subnode of {s} in {D}.

  In particular, the null state {(0,NULL)}is a substate of every
  non-accepting state, and the unit state {(1,NULL)} is a substate of
  every accepting state.

  LEXICOGRAPHIC ORDER: The ordering of the input symbols induces a
  /lexicographic ordering/ of any finite set of input strings, and hence
  of all the paths out of a state {u}.

  IMMUTABILITY: Except for {raut_crunch}, the procedures and
  methods in this interface never modify or destroy existing states,
  but merely create new ones as needed. Therefore, any properties of
  any state {u} that depend only on {u} and its successors (such as
  the suffix set) will remain constant while {raut_crunch} is not
  used.

  Moreover, {raut_crunch}, {raut_root_get_set}, {rdag_add_string},
  {rdag_remove_string}, and {rdag_build} are the only procedures in
  this interface that change the current root state. It follows the
  prefix set {Pref(u)} of any state {u} is constant as long as these
  methods are not used. */

  /* Procedures that deal with prefixes ("raut_prefs_count", {raut_enum_prefs}, ...)  and
    those that deal with incoming edges ("raut_in_degree", {raut_enum_in_arcs}) automatically
    construct internal tables whose cost is proportional to {raut_max_state(D)}.
    The tables are reused by subsequent calls to those methods, but are
    invalidated by {raut_crunch} and {raut_root_get_set}.  So, the first prefix-related
    method call after a {raut_crunch} or {raut_root_get_set} may be a lot more expensive
    than typical calls. */

/*
  DOCUMENTATION STRING: An automaton also has a documentation string,
  which may be NULL or empty, and may contain embedded line breaks
  (LF = '\n'). This string is written by {raut_dump} and
  read by {raut_load}, but is not otherwise used. */
  
char *raut_doc_get(raut_t *A);
  /* Returns the documentation string of {A}. */

void raut_doc_set(raut_t *A, char *doc);
  /* Sets the documentation string of {A} to {doc}. The previous
    string is *not* reclaimed. */
 
/*
  SPACE MANAGEMENT: */
 
rdag_node_t raut_node_max_alloc(raut_t *A);
  /* Max node for which storage has been allocated. Clients can create
    at least {raut_node_max_alloc() - raut_node_max()} new nodes
    without triggering an automatic expansion. */
 
uint32_t raut_reachable_node_count(raut_t *A, raut_state_t u);
  /* Count of all distinct proper dag nodes that are reachable from 
    state {u} by any combination of taking substates or following outgoing arcs.  The accept bit of 
    {u} is immaterial. The count includes {u.nd} itself, if not null; 
    and is zero if and only if {u.nd} is null. Cost: {O(u.nd)}
    time and space. */
 
void raut_crunch(raut_t *A);
  /* Performs an {rdag_crunch} on the dag of {A}.
    The root of {A} is updated appropriately.

    WARNING: By definition, {raut_crunch} generally changes the
    nodes of all states, and deletes some of them. */

void raut_expand(raut_t *A, rdag_node_t max_alloc_node);
  /* Expands the internal storage to accomodate nodes
    from 0 to {raut_node_max_alloc}. Preserves all current nodes
    and their numbers. Bombs out if there is no more space,
    or if {raut_node_max_alloc} exceeds the maxmum possible 
    node number. */

void raut_free(raut_t *A);
  /* Frees all storage used by the automaton {A}, except the documentation string
    but including its underlying dag and all internal tables. */

//  raut_state_t raut_max_state(raut_t *A);
//        /*
//          The number of the highest-numbered state in the automaton.
//          Note: not all numbers in {[0..raut_max_state()]} correspond
//          to currently reachable states.
//          Cost: {O(1)} time, 0 space. */
//  
//        Last(raut_t *A, raut_state_t u): Arc;
//        /*
//          The last proper arc out of {s}, in intrinsic order.
//          Requires {raut_state_has_arcs(s)}.
//          Cost: {O(1)} time, 0 space. */
//  
//        Rest(raut_t *A, raut_state_t u): State;
//        /*
//          The state obtained from {s} by removing the transition {Last(s)},
//          with the same accepting bit as {s}. Requires {raut_state_has_arcs(s)}.
//          Cost: {O(1)} time, 0 space. */
//  
//        /******************************************************************/
//        /* DERIVED METHODS                                                */
//        /******************************************************************/
//  
//        /*
//          These methods could be ordinary procedures, since
//          they can be defined in terms of the primitives above;
//          they are declared here as methods to reduce client confusion. (?) */
//  
//        OutDeg(raut_t *A, raut_state_t u): NAT;
//        /*
//          Number of proper arcs out of state {s} (0 iff {s == rdag_state_NULL}
//          or {s == UnitState}). Cost: {O("OutDeg(s)")} time, 0 space. */
//  
//        raut_in_degree(raut_t *A, raut_state_t u): NAT;
//        /*
//          Number of proper arcs that lead to state {s} from some state reachable
//          from {Root()}.  Returns 0 is {s == Root()} or {s == rdag_state_NULL}.
//          Cost: {O("raut_in_degree(s)")} time, 0 space (plus the cost of
//          rebuilding the reverse-data tables, if {Root()} has changed). */
//  
//        First(raut_t *A, raut_state_t u): Arc;
//        /*
//          The first arc out of state {s}, in intrinsic order. 
//          Requires {s!=rdag_state_NULL}. Cost: {O("OutDeg(s)")} time, 0 space. */
//  
//        Step(raut_t *A, raut_state_t u; symbol: Symbol): State;
//        /*
//          The successor of the given state through the arc labeled with
//          the given symbol.  Returns {rdag_state_NULL} if {s} is {rdag_state_NULL} or has
//          no outgoing edge labeled with that {symbol}.
//          Cost: {O("OutDeg(s)")} time, 0 space. */
//  
//        SetArc(raut_t *A, raut_state_t u; symbol: Symbol; dest: State): State RAISES {Full};
//        /*
//          The unique state that has the same arcs and accepting bit as {s},
//          except that it goes to {dest} through the given symbol.
//          Creates the state if it doesn't exist yet.
//  
//          In particular, returns {s} if {Step(s, symbol) == dest}
//          already.  The {dest} state may be {rdag_state_NULL}, in which case
//          the arc from {s} through the given symbol is removed.
//          Cost: {O("OutDeg(s)")} time and space, if hashing works. */
//  
//        SetFinal(raut_t *A, raut_state_t u; accepting: BOOLEAN): State RAISES {Full};
//        /*
//          The unique state that has the same arcs as {s}, except that
//          its accepting bit is {accepting}.  Creates the state if it doesn't
//          exist yet.
//  
//          The resulting state accepts all the strings accepted by {s},
//          plus the empty string.  In particular, the call
//          {SetFinal(rdag_state_NULL, TRUE)} returns unit state.
//  
//          Cost: {O("OutDeg(s)")} time and space, if hashing works.
//          Hence, when building a state from scratch, it is better to set
//          the accepting bit before adding the outgoing edges. */
//  
//        Walk(raut_t *A, raut_state_t u; READONLY w: String): State;
//        /*
//          Spells {w} in the automaton starting from {s}, that is,
//          returns the state reached from {s} by the unique path
//          whose edges are labeled with the symbols of {w}. */
//  
//        Accepts(raut_t *A, raut_state_t u; READONLY w: String): BOOL;
//        /*
//          TRUE iff {w} is in the language {Suff(s)}.  Equivalent to
//          {Accepting(Walk(s, w))}. */
//  
//        Rank(raut_t *A, raut_state_t u; READONLY w: String; reverse: BOOL = FALSE): NAT;
//        /*
//          The rank of {w} in {Suff(s)}, in lexicographic order.  More
//          precisely, the number of words in {Suff(s)} that are strictly
//          smaller than {w} (or strictly greater, if {reverse == TRUE}).
//          The result is defined even if {w} is not in {Suff(s)}.  Cost:
//          {O*(n(d + 1))}, where {n} is the length of {w}, and {d} is the
//          average value of {OutDeg(t)} over all states {t} in the path
//          spelled by {w}. */
//  
//        AddString(raut_t *A, raut_state_t u; READONLY w: String): State RAISES {Full};
//        /*
//          Returns the state {s'} such that {Suff(s') == Suff(s) + {w}}.
//          (In particular, returns {s} itself if {w} is already in {Suff(s)}.
//  
//          Note that {AddString} does NOT modify the current root
//          state, and does not change {Pref(t)} for any state {t}. */
//  
//        SubString(raut_t *A, raut_state_t u; READONLY w: String): State RAISES {Full};
//        /*
//          Returns the state {s'} such that {Suff(s') == Suff(s) - {w}}.
//          In particular, returns {s} itself if {w} is not in Suff(s).
//  
//          Note that SubString does NOT modify the current root
//          state, and does not change {Pref(t)} for any state {t}.
//  
//          In the current implementation, the spacetime cost of
//          {AddString} and {SubString} is {\Theta(m n)} in the worst case,
//          where {m} is the number of symbols in the automaton's alphabet,
//          and {n} is the length of {w}. */
//          
//        AddStringMayberaut_crunch(raut_t *A, 
//            State *s; 
//            String READONLY w; 
//            keep: REF States = NULL;
//          ): State;
//        
//        SubStringMayberaut_crunch(raut_t *A, 
//            State *s;
//            String READONLY w; 
//            keep: REF States = NULL;
//          ): State;
//        /*
//          Like {AddString} and {SubString}, except that if there isn't enough
//          space in the automaton, calls {raut_crunch(s & keep)}, and tries again.
//          If {raut_crunch} doesn't return enough space, expands the automaton
//          and tries again, until it succeeds (or runs out of memory). */
//          
//        /***************************************************************************/
//        /* ENUMERATION                                                             */
//        /***************************************************************************/
//  
//        EnumOutArcs(raut_t *A, raut_state_t u; action: ArcAction) RAISES {Skip, Abort};
//        /*
//          Enumerates the arcs out of state {s}, in the intrinsic order
//          (from {First(s)} to {Last(s)}), and calls {action(i,a)} for each;
//          where {i} is the index of the arc (counting from 0) and {a} the
//          Arc data. 
//  
//          If {action(i, a)} raises {Skip}, {EnumOutArcs} ignores all
//          arcs out of {s} that follow {a} in the intrinsic ordering.
//          {EnumOutArcs} passes to the client any {Abort}s raised by
//          {action}.*/
//  
//        raut_enum_in_arcs(raut_t *A, raut_state_t u; action: ArcAction) RAISES {Skip, Abort};
//        /*
//          Enumerates the arcs that enter {s}, and calls {action(i, a)} for
//          each arc, where {i} is the index of the arc in the enumeration
//          (counting from 0), and {a} is the data for the arc with its
//          direction reversed (i.e., {a.dest} is actually the origin of 
//          the arc; the actual destination is always {s}).
//  
//          The arcs are visited in order of increasing origin state
//          {a.dest}; arcs with same origin {t} are visited in intrinsic 
//          order, from {First(t)} to {Last(t)}.
//  
//          The exceptions {Skip} and {Abort} are handled just as in {EnumOutArcs}. */
//  
//        EnumPaths(raut_t *A, 
//            State s;
//            enter: StateAction = NULL;
//            push: PathAction = NULL;
//            pop: PathAction = NULL;
//            exit: StateAction = NULL;
//          ) RAISES {Abort};
//        /*
//          Enumerates paths starting from {s}, in standard double order.
//          Only paths that consist entirely of proper arcs and states are
//          considered.
//          
//          The method actually enumerates only a subset of all such
//          paths, as determined by the client-provided procedures {enter,
//          }push", {pop}, and {exit}.
//  
//          The enumeration algorithm can be described in terms of a conceptual
//          {current path} that initially contains just the state {s},
//          and grows or shrinks one arc at a time.  {EnumPaths} calls
//  
//              {enter(len,o,f)} when the path has {len} arcs and has just landed
//                     on a state {o}, whose {accepting} bit is {f}.
//  
//              {push(len,o,i,a)} whenever the current path has {len} arcs,
//                     ends at state {o}, and is being extended with the
//                     arc {a}, which is the {i}th arc out of {o} (where
//                     {i == 0} means arc {First(o)}).
//  
//              {pop(len,o,i,a)} whenever the current path has {len+1} arcs
//                     and EnumPaths is about to remove its last arc {a},
//                     which is the {i}th arc out of {o}.
//  
//              {exit(len,o,f)} whenever the current path has {len} arcs and
//                     ends in state {o} with accepting bit {f}, and no extensions
//                     of it remain to be enumerated.
//  
//          Any of the four client actions can stop the enumeration by raising {Abort}.
//  
//          Unless the enumeration is aborted, every call to {enter} will be followed
//          eventually by a matching call to {exit}, and every call to {push}
//          will be followed eventually by a matching call to {pop}.  The typical
//          action pattern for a generic node {t} with {n} outgoing proper arcs is
//  
//          |     enter(len,t,f)
//          |       push(len,t,0,a)
//          |         ...
//          |       pop(len,t,0,a)
//          |       push(len,t,1,b)
//          |         ...
//          |       pop(len,t,1,b)
//          |       ...
//          |       push(len,t,n-1,x)
//          |         ...
//          |       pop(len,t,n-1,x)
//          |     exit(len,t,f)
//  
//          In particular, if the starting state {s} is non-null, {EnumPaths} will call
//          {enter(0,s,Accepting(s))} at the very beginning, and (if not aborted)
//          will call {exit(0,s,Accepting(s))} at the very end.
//          (If the starting state {s} is {rdag_state_NULL}, {EnumPaths} will do nothing.)
//  
//          Note that the same state {o} may be {enter}ed and {exit}ed many, many times.
//  
//          The client can prune branches of the path tree by raising the
//          {Skip} exception in the action procedures:
//  
//              if {enter(len,o,f)} raises {Skip}, {EnumPaths} will omit the 
//              entire branch of the path tree rooted at current path, and call 
//              {exit(len,o,f)} right away, as if {o} had no outgoing arcs;
//  
//              if {push(len,o,i,a)} raises {Skip}, {EnumPaths} will
//              call {pop(len,o,i,a)} right away, as if the arc {a} led
//              to {rdag_state_NULL};
//  
//              if {pop(len,o,i,a)} raises {Skip}, {EnumPaths} ignores any remaining
//              arcs out of {o}, and calls {exit(len,o,f)} right away.
//  
//          (It is a checked error for {exit(len,o,f)} to raise {Skip}.)
//  
//          If any action is not specified, {EnumPaths} will provide a
//          a trivial default action that does nothing (and raises no exceptions).
//          */
//  
//        EnumStrings(raut_t *A, raut_state_t u; enter, exit: StringAction = NULL) RAISES {Abort};
//        /*
//          Enumerates strings spelled by paths out of {s}, in standard double
//          order.
//  
//          {EnumStrings} is similar {EnumPaths}, but gives to the {enter} and {exit}
//          actions the whole string spelled by the current path, and its last
//          state.
//  
//          Specifically, {EnumStrings} calls
//  
//              {enter(w, d, f)} when it has reached for the first time
//                  some state {d}, whose {accepting} bit is {f}, by spelling string {w},
//                  and it is about to enumerate all suffixes of {s} whose 
//                  proper prefix is {w};
//  
//              {exit(w, d, f)} when it is going to leave state {d}, whose
//                  {accepting} bit {f}, and has just finished enumerating
//                  all suffixes of {s} whose proper prefix is {w}.
//  
//          Either {enter} or {exit} can stop the enumeration by raising {Abort}.
//  
//          Unless the enumeration is aborted, every call to {enter} will be followed
//          eventually by a matching call to {exit}.  The very first call to {enter}
//          will be {enter({}, s, f)}, and the very last call to {exit} will
//          be {exit({}, s, f)}, where {{}} is the empty input strings.
//  
//          If {enter(w, d, f)} raises {Skip}, {EnumStrings} will omit the entire branch
//          of the path tree rooted at the resulting current path, and will call
//          {exit(w, d, f)} right away.   (It is an error for {exit} to raise {Skip}.)
//  
//          Note that {EnumStrings} enumerates all the strings that lead from {s}
//          to *any* proper state, not just those that lead to a accepting state.
//          If you only want the latter, use {EnumSuffs}.
//  
//          If either {enter} or {exit} is not specified,  {EnumStrings} supplies
//          a trivial default {StringAction} that does nothing (raut_t *A, and raises no exceptions).
//          */
//  
//        EnumStates(raut_t *A, 
//            States READONLY base;
//            substates: BOOL = FALSE;
//            enter: StateAction = NULL;
//            exit: StateAction = NULL;
//          ) RAISES {Abort};
//        /*
//          The procedure enumerates all states reachable from the states
//          base[i], and calls the client procedures {enter} and {exit}
//          exactly once each for every non-null state it sees.
//  
//          If {substates == FALSE} (the default), the procedure only
//          enumerates those states that can be reached from {base} by
//          zero or more {Step} steps. If {substates == TRUE}, the procedure
//          also enumerates all substates of those states; in other words,
//          all states that can be reached from {base} by
//          zero or more {Step} or {Rest} operations.
//          In either case, {EnumStates} visits each state exactly once,
//          even if it is a substate shared by two or more reachable states.
//  
//          The states are enumerated twice in numerical order (which is consistent
//          with both the {Step} and {Rest} operations): once in decreasing
//          order with {enter}, and once in increasing order with {exit}.
//  
//          More precisely, {EnumStates} will call {enter(len,s,f)} for each reachable
//          state {s}, in *decreasing* numerical order, where {len} is the minimum
//          number of {Step} steps (or {Step} and {Rest} steps, if {substates==TRUE})
//          from {base} to {s}, and {f} is {Accepting(s)}.  Note that this call will happen
//          before any successors or substates of {s} are {enter}ed.
//  
//          That done, {EnumStates} will call {exit(len,s,f)} for each reachable
//          state {s}, in *increasing* numerical order, where {len} and {f} have the
//          same meaning as before.  Note that this call will happen after
//          all successors and reachable substates of {s} have been {exit}ed.
//  
//          The actions {enter} and {exit} can stop the enumeration
//          at any time by raising {Abort}.
//  
//          The {enter} procedure can also control the set of visited vertices through
//          its returned value. If {enter(len,s,f)} raises {Skip}, {EnumStates} will
//          ignore any arcs out of {s}.  In that case, the states that can be reached
//          from {s} will be visited only if they are reachable through other paths
//          that do not include {s}.  (It is an error for {exit} to raise {Skip}.)
//  
//          If either {enter} or {exit} is not specified, {EnumStates} will
//          supply a trivial {StateAction} that does nothing (and raise no exceptions)
//     .
//          */
//  
//        /***************************************************************************/
//        /* PREFIXES)  AND  AND  (SUFFIXES                                                   */
//        /***************************************************************************/
//  
//        raut_prefs_count(raut_t *A, raut_state_t u): NAT;
//        /*
//          The number of distinct paths that lead from the current
//          root state to the state {s}. The result is zero
//          if the current root is {rdag_state_NULL}, if {s} is rdag_state_NULL,
//          or {s} is not reachable from the current root. */
//  
//        NSuffs(raut_t *A, raut_state_t u): NAT;
//        /*
//          The number of distinct paths leading from state {s} to some
//          accepting state. The result is zero iff {s} is {rdag_state_NULL}. */
//  
//        NPrefLetters(raut_t *A, raut_state_t u): NAT;
//        /*
//          The total length of all prefixes of {s}, that is,
//          of all paths leading from the current root state to state {s}.
//          The result is zero if the current root is {rdag_state_NULL}, if {s} 
//          is {rdag_state_NULL}, or {s} is not reachable from the current root. */
//  
//        NSuffLetters(raut_t *A, raut_state_t u): NAT;
//        /*
//          The total length of all suffixes of {s}, that is,
//          of all paths leading from {s} to some accepting state.
//          The result is zero iff {s} == rdag_state_NULL or {s == UnitState}. */
//  
//        FirstPrefix(raut_t *A, raut_state_t u; e: Encoding.T): char *;
//        /*
//          The first string (in the funny order defined by {raut_enum_prefs})
//          leading from the root to {s}, converted to a char *by applying
//          {lpr} to each symbol.   An error if {s} is {rdag_state_NULL}, or
//          is not reachable from the current root (in particular,
//          if the current root is {rdag_state_NULL}). */
//  
//        FirstSuffix(raut_t *A, raut_state_t u; e: Encoding.T): char *;
//        /*
//          The first string (in standard pre-order) leading from {s}
//          to some accepting state, converted to char *by applying {lpr} to
//          each symbol.  An error if {s} is {rdag_state_NULL}. */
//  
//        FullLabel(raut_t *A, raut_state_t u; e: Encoding.T; sep: char *= {:}): char *;
//        /*
//          Computes a label for the state {s}, of the form
//          {<first-prefix><sep><first-suffix>}.
//          An error if {s} is not reachable from the root,
//          or if {s} is {rdag_state_NULL}. */
//  
//        raut_enum_prefs(raut_t *A, raut_state_t u; action: PrefixAction) RAISES {Abort};
//        /*
//          Enumerates {Pref(s)}, the set of all strings leading from the
//          root to state {s}, calling {action} on each.
//          The state {s} must not be {rdag_state_NULL}.
//  
//          IMPORTANT (1): each string passed to {action} is actually
//          the symbol-reversal of the prefix in question.
//  
//          IMPORTANT (2): The strings are enumerated in a funny order
//          (depth-first order based on {raut_enum_in_arcs}) that is NOT
//          lexicographic order, and is not related to the standard
//          path orders.
//  
//          {raut_enum_prefs} will return without calling {action}
//          if {s} is not reachable from the current root.
//  
//          The {action} may raise {Abort} to stop the enumeration. */
//  
//        EnumSuffs(raut_t *A, raut_state_t u; action: SuffixAction) RAISES {Abort};
//        /*
//          Enumerates {Suff(s)}, the set of all strings leading from
//          {s} to some accepting state, calling {action} on each string.
//  
//          The strings are enumerated in standard path pre-order.
//  
//          {EnumSuffs} will return without calling {action}
//          iff {s} is {rdag_state_NULL}.
//  
//          The {action} may raise {Abort} to stop the enumeration. */
//  
//        PrintPrefs(raut_t *A, raut_state_t u; spr: StringPrinter.T);
//        /*
//          Prints all prefixes of state {s}, in the funny order used by {raut_enum_prefs}.
//          Uses the given {StringPrinter.T} for formatting control.
//          An error if {s} is {rdag_state_NULL} or is not reachable from the root. */
//  
//        PrintSuffs(raut_t *A, raut_state_t u; spr: StringPrinter.T);
//        /*
//          Prints the suffixes of state {s}, in standard path pre-order.
//          Uses the given {StringPrinter.T} for formatting control.
//          Prints nothing iff {s} is {rdag_state_NULL}. */
//  
//        /***************************************************************************/
//        /* CONSTRUCTION                                                            */
//        /***************************************************************************/
//  
//        rdag_build(raut_t *A, 
//            NextStringProc next;        /* Client input procedure */
//            wr: Wr.T = NULL;             /* Writer for messages, progress report, etc. */
//            e: Encoding.T = NULL;        /* Symbol/CHAR encoding for messages */
//            reportInterval: POS = 1000; /* Status reporting interval */
//            flagRedundant: BOOL = TRUE; /* TRUE to print warnings on redundant ops */
//          ) RAISES {Abort};
//        /*
//          Modifies the automaton by adding and/or deleting given set of strings.
//          
//          {rdag_build} calls {next(s, len, add)} repeatedly to get each string.  
//          The {next} procedure should return the string's length in {len}, and
//          the string itself in {s^[0..len-1]} (expanding s^ if necessary).  
//          Also, {next} should set {add = TRUE} for strings to be
//          added, and {add = FALSE} for strings to be deleted.  Finally, {next}
//          should raise {Done} when there is no next string.
//  
//          {rdag_build} will automatically crunch and/or expand the automaton as necessary.
//          
//          If {wr} is not NULL, {rdag_build} will periodically print to {wr} a one-line
//          progress report, giving the string and symbol count and the size of
//          the automaton.  The report will be printed after processing
//          {reportInterval} additional strings, after each {raut_crunch} and {Expand},
//          and before returning. 
//          
//          The Encoding {e} is used only to print out strings in error messages.
//          It defaults to a {PlainEncoding.New()}. */
//  
//        /***************************************************************************/
//        /* STATISTICS                                                              */
//        /***************************************************************************/
//  
//        NStates(raut_t *A, READONLY base: States): NAT;
//        /*
//          The number of distinct states reachanble from {base} by {Step} paths. 
//          Cost: {O(N)}. */
//  
//        NArcs(raut_t *A, READONLY base: States): NAT;
//        /*
//          The number of distinct arcs reachable from the {base} states by 
//          {Step} paths.  Cost: O("raut_max_state()"). */
//  
//        Count(raut_t *A, READONLY base: States): Counts;
//        /*
//          Computes:
//  
//            o the number of distinct non-null states, accepting states, and arcs 
//              in the part of {self} that is reachable from the given {base} 
//              states by {Step} paths.
//  
//            o the number of distinct non-null substates reachable from {base}
//              by {Step}+"Rest" paths;
//  
//            o the total count and length of all strings accepted by 
//              the states in {base}.
//  
//          Every reachable state, arc, accepting, or sub-state is counted
//          only once, irrespective of how many {base} elements lead to it.
//          On the other hand, if two {base} states accept the same string
//          (in particular, if they are the same state), that string
//          contributes twice to the string and length counts.
//          */
//  
//        /***************************************************************************/
//        /* STORAGE MANAGEMENT                                                      */
//        /***************************************************************************/
//  
//        MaxAllocState(raut_t *A, ): POS;
//        /*
//          Implementation data: maximum number of states/substates
//          for which storage has been allocated.
//  
//          In the current implementation, each distinct non-null state
//          or substate  (reachable or not) occupies one unit of storage.
//          Attempts to create more than {MaxAllocState()} distinct substates
//          will generally raise the {Full} exception. */
//  
//        Expand(raut_t *A, newSize: NAT = 0) RAISES {Full};
//        /*
//          Reallocates enough storage for up to {newSize} new substates.
//          Preserves all current states and their numbers.
//          If {newSize} is not given, uses a constant factor
//          (between 1 and 2) times the current size. Raises {Full}
//          iff the automaton already contains more than {newSize}
//          distinct non-null substates. */
//  
//        raut_crunch(raut_t *A, keep: REF States = NULL);
//        /*
//          Discards all states not reachable from the current root
//          or from the states listed in {keep} (if any); then
//          squeezes the reachable ones together as closely as possible,
//          and updates the current root and the {keep} vector to reflect
//          the new State numbering.
//  
//          WARNING: By definition, {raut_crunch} generally changes the numbers
//          of all states, and deletes some of them.  The client should
//          copy all its important State variables into the {keep} vector,
//          and copy them back after the {raut_crunch} is done. */
//  
//    ;};
//  
//  TYPE
//    ArcAction == PROCEDURE(i: NAT; a: Arc) RAISES {Skip, Abort};
//    PathAction == PROCEDURE(len: NAT; org: State; i: NAT; a: Arc) RAISES {Skip, Abort};
//    StateAction == PROCEDURE(len: NAT; raut_state_t u; accepting: BOOL) RAISES {Skip, Abort};
//    StringAction == PROCEDURE(READONLY w: String; raut_state_t u; accepting: BOOL) RAISES {Skip,Abort};
//    PrefixAction == PROCEDURE(READONLY w: String) RAISES {Abort};
//    SuffixAction == PROCEDURE(READONLY w: String) RAISES {Abort};
//  
//    NextStringProc == PROCEDURE(
//      VAR /*IO*/ s: REF String;  /* The string buffer */
//      VAR /*OUT*/ len: NAT;      /* Length of returned string */
//      VAR /*OUT*/ add: BOOLEAN;  /* TRUE to add, FALSE to delete */
//    ) RAISES {Done, Abort};
//  
//  TYPE
//    Counts == RECORD
//        NAT strings;   /* Total strings accepted from all roots */
//        NAT symbols;   /* Total length os all strings accepted from all roots */
//        NAT states;    /* Distinct non-null states reachable by {Step} */
//        NAT finals;    /* Distinct reachable states that are accepting */
//        NAT arcs;      /* Total number of outgoing arcs from all reachable states */
//        NAT substates; /* Distinct non-null substates reachable by {Step}+"Rest" */
//      ;};
//  
//  /***************************************************************************/
//  /* CREATION                                                                */
//  /***************************************************************************/
//  
//  PROCEDURE Copy(
//      T aut;
//      size: POS = 1;
//      VAR /*OUT*/ map: REF States;
//    ): T RAISES {Full};
//  /*
//    Returns a new automaton that contains {rdag_state_NULL}, unit state, and an isomorphic
//    copy of the part of {aut} that is reachable from {aut.Root()} by chains of
//    {Step} or {Rest} (i.e. all reachable states of {aut}, and their substates).
//  
//    The procedure also returns a vector {map} that gives the correspondence
//    between old and new state numbers: specifically, {map[s]} is the state of the
//    copy that corresponds to state {s} of {aut}.  If {s} has no corresponding
//    state in the copy, the preocedure sets {map[s]} to rdag_state_NULL.
//  
//    The {size} parameter has the same meaning as in {New}.
//    The procedure may raise {Full} if {size} is not enough.
//    */
//  
//  PROCEDURE CopyStates(
//      T from;
//      T to;
//      State s;
//      VAR /*IO*/ map: REF States;
//    ): State RAISES {Full};
//  /*
//    Copies into the {to} automaton all states of the {from} automaton
//    that are reachable from {root} by {Step} chains.
//  
//    If the {map} argument is not NULL , {Copy} assumes that any state {t} with
//    {map[t]!=rdag_state_NULL} has already been copied into the {to} automaton, and its
//    number there is {map[t]}.  In that case, {Copy} will set {map[t]}
//    appropriately for every state that it copies.  (Note that if {map} is
//    non-NULL, one must have {LAST(map^) >= s}. */
//  
//  /***************************************************************************/
//  /* IO                                                                      */
//  /***************************************************************************/
//  
//  PROCEDURE Dump(wr: Wr.T; aut: T);
//  PROCEDURE Load(rd: Rd.T; minSize: POS = 1): T;
//  /*
//    {Dump} writes {aut} to the given writer, in a format that can be
//    read back with {Load}.  
//    
//    The output is in ASCII, but it is not meant to be human-readable. (except for
//    the {doc} text, which is written at the beginning of the file.)
//  
//    {Dump} only writes the existing states and arcs, from {1} to
//    {aut.raut_max_state()}. The {Reduced.T} returned by {Load} is just large enough to
//    contain those arcs, or {minSize} arcs, whichever is larger. */
//  
//  PROCEDURE Print(wr: Wr.T; aut: T; e: Encoding.T);
//  /*
//    Writes to {wr} a legible description of the automaton states that
//    are reachable from the current root state, in the format
//  
//      | root == nnnn
//      | 
//      | 418 *
//      |   a -> 4615
//      |   b -> 4634
//      |   h -> 1950
//      | 
//      | 4615
//      |   a -> 5019
//      |   ...
//  
//    Each paragraph describes one state. The first line of each paragraph gives
//    the state number, and a '*' iff the state is accepting.  The remaining lines give
//    the symbol labels and the destinations of the arcs out of the state.
//  
//    The input symbols label of each arc is printed using the procedure
//    {e.PrintLetter}. */
//  
//  PROCEDURE PrintCounts(wr: Wr.T; READONLY ct: Counts);
//  /*
//    Prints the given counts to {wr}. */

#endif
