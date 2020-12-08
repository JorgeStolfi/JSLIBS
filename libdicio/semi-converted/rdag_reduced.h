#ifndef acau_reduced_H
#define acau_reduced_H

/* acau_last edited on 2009-10-23 21:48:22 by stolfi */

#define acau_reduced_DESC "Representation of reduced deterministic acyclic automata"

/*
  An {acau_reduced_t} is a deterministic acyclic finite automaton,
  reduced but not necessarily connected.

  STATES, ARCS: an {acau_reduced_t} is a directed graph with a finite number
  of vertices, the /proper states/.  Each state has zero or more
  outgoing directed edges, the /proper arcs/ of that state,
  each labeled with a distinct symbol and leading to another proper state.

  FINAL BIT: each state has a {final} attribute that may be TRUE or FALSE.

  ACYCLICITY: the directed graph defined by the proper states and proper
  arcs of a {acau_reduced_t} is acyclic, and from every state there is at least
  one path leading to a final state.

  NULL STATE: it is convenient to assume that every {acau_reduced_t} has,
  in addition to the proper states described above, a {virtual} state
  {acau_red_state_NULL} with {final==FALSE} and no arcs.  It is also convenient
  to assume that whenever a state has no proper transition labeled
  with a given symbol, it has a /virtual/ transition labeled with
  that symbol and leading to {acau_red_state_NULL}.  In particular, the
  {acau_red_state_NULL} has a virtual transition to itself under any symbol.
  (Thus, the graph is NOT acyclic if {acau_red_state_NULL} and it virtual arcs
  are included.)

  ROOT STATE: A {acau_reduced_t} has a distinguished /root state/
  {acau_red_root()}, which affects the meaning of certain operations
  (notably those having to do with prefix strings, see below).

  In general, there may be states in the automaton that are not reachable
  from its current root.  In fact, there may be no single state from
  which all other states are reachable.  In particular, the root state
  may be {acau_red_state_NULL}.

  SUFFIXES OF A STATE: each proper state {s} is said to {accept} or {recognize}
  a language {acau_red_suff(s)}, the set of /suffixes/ of {s}, which is the set of all
  {String}s spelled by directed paths that lead from {s} to some final state.
  Obviously, the the language {acau_red_suff(s)} for a proper state is always finite
  and non-empty.  By definition, {acau_red_suff(acau_red_state_NULL)} is the empty set.

  PREFIXES OF A STATE: each proper state {s} has also a set of /prefixes/
  {acau_red_pref(s)}, the set of all {String}s spelled by paths that lead from the
  root state to {s}.  In particular, if {s} is the current root state,
  then {acau_red_pref(s)} contains only the empty string; and if {s} is not reachable
  from the root state, then {acau_red_pref(s)} is the empty set.
  
  For technical and logical reasons, the set {acau_red_pref(acau_red_state_NULL)} is undefined.

  STATE UNIQUENESS: the implementation ensures that the automaton
  is /reduced/, that is, no two distinct states (proper or null)
  have the same set of suffix strings.  In particular, the {acau_red_state_NULL} is
  the only state that accepts the empty set of strings {}.

  DETERMINISM: according to the definitions above, the outgoung arcs
  of any given state are labeled with distinct symbols.  It follows
  that the automaton is deterministic: given any state {s} and any
  String {w}, there is exactly one state {t} (possibly null) that
  can be reached from {s} by a path whole labels spell {w}.
  Therefore, any two distinct states have disjoint prefix sets;
  and any proper reachable state can be uniquely identified by any
  of its prefix strings.

  UNIT STATE: besides the {acau_red_state_NULL}, every {acau_reduced_t} has
  exactly one proper state, the {acau_red_state_UNIT}, that has {final==TRUE} and no
  outgoing arcs.  Obviously, that state accepts the language consisting
  only of the empty string {()}.
  
  INTRINSIC ARC ORDER: for every state {s} that has outgoing arcs
  (that is, any state other than {acau_red_state_NULL} or {acau_red_state_UNIT}), there is
  one distinguished arc {acau_last(s)} out of {s}, and another state
  {acau_rest(s)}, such that the set of arcs leaving the latter are those
  leaving {s}, minus the arc {acau_last(s)}.
  
  The methods {acau_last} and {acau_rest} define the /intrinsic ordering/ for the
  arcs out {s}.  Note that {acau_rest} moves {\em backwards} in this ordering.
    
  SUBSTATES: A {substate} of a state {s} of a acau_reduced_t is a state {t}
  that whose transitions are an initial prefix of the transitions out
  of {s}, in intrinsic order (or, in other words, that can be obtained
  from {s} by zero or more applications of the {acau_rest} method).

  In particular, the {acau_red_state_NULL} is a substate of every state
  with {final==FALSE}, and the {acau_red_state_UNIT} is a substate of every
  state with {final==TRUE}.

  If a {acau_reduced_t} includes a state {s}, it automatically includes
  all substates of {s}.
  
  TOPOLOGICAL ORDERING OF THE STATES: The states in a {acau_reduced_t} are
  identified by numbers, assigned as they are created.  The numerical
  order of the states is compatible with the transitions and the
  substate relation; that is, every state is greater (numerically)
  than any of its proper substates, and also greater than the
  destination of any of its outgoing arcs.

  STANDARD PATH ORDERS: The intrinsic ordering of the arcs out of each
  state induces a /standard/ order for the set of all paths out of a
  fixed state {s}.  In this order, all paths that leave {s} through
  the arc {acau_last(s)} occur after those that leave through other arcs.
  
  There are actually three standard path orderings, differing on where
  the paths that {\em end} at {s} appear: either before all the paths
  that {\em go through} {s} (/pre-order/), or after them
  (/end-order/), or both before and after them (/double order/).
  
  STRING SPELLING: Every path out of a state {s} is said to
  /spell out/ the {String} consisting of the symbols that label the
  arcs along the path.  Conversely, since the automaton is
  deterministic, each {String} defines a unique path out of a given
  state {s}.
  
  LEXICOGRAPHIC ORDER: The ordering of the {acau_symbol_t} values induces a
  /lexicographic ordering/ of any finite set of {String}s, and hence
  of all the paths out of a state {s}
  
  Note that the lexicographic ordering is distinct from all the
  standard path orderings (unless the arcs in every state {s} happen
  to be sorted by symbol value, with {acau_last(s)} being the highest).

  IMMUTABILITY:  Except for the {acau_red_crunch} method (below), the procedures
  and methods in this interface never modify or destroy existing states,
  but merely create new ones as needed.  Therefore, if {acau_red_crunch} is not used,
  the {final} bit and outgoing arcs of a state {s} are constant
  throughout its lifetime, and so are any properties of {s}
  that can defined in terms of these (such as the suffix set).

  Moreover, {acau_red_crunch}, {acau_red_set_root}, and {acau_red_build} are the only two methods or
  procedures in this interface that change the current root state.  It follows
  the prefix set acau_red_pref(s) of any state {s} is constant as long as these methods
  are not used.  */

#include <acau_basics.h, Rd, Wr, acau_string_printer.h, acau_encoding.h>
#include <acau_basics.h> /* uint32_t, POS, BOOL; */
#include <acau_basics.h> /* Done, acau_red_error_FULL, acau_disp_SKIP, acau_disp_STOP; */

/*
  The methods and procedures below raise {acau_basics.h.acau_red_error_FULL} if there isn't
  enough allocated storage in the automaton for the requested operation. */

typedef
  acau_symbol_t == [ FIRST(acau_basics.h.acau_symbol_t)+1 .. LAST(acau_basics.h.acau_symbol_t) ];
  /*
    Note that FIRST(acau_basics.h.acau_symbol_t) is reserved for internal
    use of the implementation */

  String == acau_basics.h.String;
  /*
    Should be ARRAY OF Reduced.acau_symbol_t, excluding FIRST(acau_basics.h.acau_symbol_t),
    but that would require special versions of ExpandString, etc. */

  acau_red_state_t == uint32_t;
  States == ARRAY OF acau_red_state_t;
  
  Arc == RECORD symbol: acau_symbol_t; dest: acau_red_state_t ;};
  acau_red_arc_vec_t == ARRAY OF Arc;

CONST
  acau_red_state_NULL == 0; /* The unique (improper) state with no arcs and {final==FALSE} */
  acau_red_state_UNIT == 1; /* The unique (proper) state with no arcs and {final==TRUE} */

typedef
  T <: Public;

  Public == OBJECT
  
      char *doc;  
        /* 
          An arbitrary multi-line text that describes the automaton. 
          If not empty, it should end with a newline character. */
      
    METHODS

      /*
        Note: the methods that deal with prefixes ("acau_red_num_prefs", {acau_red_enum_prefs}, ...)  and
        those that deal with incoming edges ("acau_red_in_deg", {acau_red_enum_in_arcs}) automatically
        construct internal tables whose cost is proportional to {acau_red_max_state()}.
        The tables are reused by subsequent calls to those methods, but are
        invalidated by {acau_red_crunch} and {acau_red_set_root}.  So, the first prefix-related
        method call after a {acau_red_crunch} or {acau_red_set_root} may be a lot more expensive
        than ordinary calls. */

      acau_red_max_state(): POS;
      /*
        The number of the highest-numbered state in the automaton.
        Note: not all numbers in {[0..acau_red_max_state()]} correspond
        to currently reachable states.
        Cost: $O(1)$ time, 0 space. */

      acau_red_is_final(s: acau_red_state_t): BOOL;
      /*
        TRUE iff {s} is a final state of the automaton.
        Cost: $O("acau_out_deg(s)")$ time, 0 space. */

      acau_red_has_arcs(s: acau_red_state_t): BOOL;
      /*
        TRUE iff the state {s} has at least one proper outgoing arc.
        Cost: $O(1)$ time, 0 space. */

      acau_last(s: acau_red_state_t): Arc;
      /*
        The last proper arc out of {s}, in intrinsic order.
        Requires {acau_red_has_arcs(s)}.
        Cost: $O(1)$ time, 0 space. */

      acau_rest(s: acau_red_state_t): acau_red_state_t;
      /*
        The state obtained from {s} by removing the transition {acau_last(s)},
        with the same final bit as {s}. Requires {acau_red_has_arcs(s)}.
        Cost: $O(1)$ time, 0 space. */

      acau_red_root(): acau_red_state_t;
      /*
        Returns the current root state.
        Cost: $O(1)$ time, 0 space. */

      acau_red_set_root(s: acau_red_state_t);
      /*
        Makes {s} the current root state.  Cost: $O(1)$ time, 0 space.
        However, note that this operation may implicitly change the
        prefix sets of many states. */

      /******************************************************************/
      /* DERIVED METHODS                                                */
      /******************************************************************/

      /*
        These methods could be ordinary procedures, since
        they can be defined in terms of the primitives above;
        they are declared here as methods to reduce client confusion. (?) */

      acau_out_deg(s: acau_red_state_t): uint32_t;
      /*
        Number of proper arcs out of state {s} (0 iff {s == acau_red_state_NULL}
        or {s == acau_red_state_UNIT}). Cost: $O("acau_out_deg(s)")$ time, 0 space. */

      acau_red_in_deg(s: acau_red_state_t): uint32_t;
      /*
        Number of proper arcs that lead to state {s} from some state reachable
        from {acau_red_root()}.  Returns 0 is {s == acau_red_root()} or {s == acau_red_state_NULL}.
        Cost: $O("acau_red_in_deg(s)")$ time, 0 space (plus the cost of
        rebuilding the reverse-data tables, if {acau_red_root()} has changed). */

      acau_first(s: acau_red_state_t): Arc;
      /*
        The first arc out of state {s}, in intrinsic order. 
        Requires {s!=acau_red_state_NULL}. Cost: $O("acau_out_deg(s)")$ time, 0 space. */

      acau_red_step(s: acau_red_state_t; symbol: acau_symbol_t): acau_red_state_t;
      /*
        The successor of the given state through the arc labeled with
        the given symbol.  Returns {acau_red_state_NULL} if {s} is {acau_red_state_NULL} or has
        no outgoing edge labeled with that {symbol}.
        Cost: $O("acau_out_deg(s)")$ time, 0 space. */

      acau_red_set_arc(s: acau_red_state_t; symbol: acau_symbol_t; dest: acau_red_state_t): acau_red_state_t RAISES {acau_red_error_FULL};
      /*
        The unique state that has the same arcs and final bit as {s},
        except that it goes to {dest} through the given symbol.
        Creates the state if it doesn't exist yet.

        In particular, returns {s} if {R(s, symbol) == dest}
        already.  The {dest} state may be {acau_red_state_NULL}, in which case
        the arc from {s} through the given symbol is removed.
        Cost: $O("acau_out_deg(s)")$ time and space, if hashing works. */

      acau_red_set_final(s: acau_red_state_t; final: BOOLEAN): acau_red_state_t RAISES {acau_red_error_FULL};
      /*
        The unique state that has the same arcs as {s}, except that
        its final bit is {final}.  Creates the state if it doesn't
        exist yet.

        The resulting state accepts all the strings accepted by {s},
        plus the empty string.  In particular, the call
        {acau_red_set_final(acau_red_state_NULL, TRUE)} returns {acau_red_state_UNIT}.

        Cost: $O("acau_out_deg(s)")$ time and space, if hashing works.
        Hence, when building a state from scratch, it is better to set
        the final bit before adding the outgoing edges. */

      acau_red_walk(s: acau_red_state_t; READONLY w: String): acau_red_state_t;
      /*
        Spells {w} in the automaton starting from {s}, that is,
        returns the state reached from {s} by the unique path
        whose edges are labeled with the symbols of {w}. */

      acau_red_accepts(s: acau_red_state_t; READONLY w: String): BOOL;
      /*
        TRUE iff {w} is in the language {acau_red_suff(s)}.  Equivalent to
        {acau_red_is_final(acau_red_walk(s, w))}. */

      acau_red_rank(s: acau_red_state_t; READONLY w: String; reverse: BOOL = FALSE): uint32_t;
      /*
        The rank of {w} in {acau_red_suff(s)}, in lexicographic order.  More
        precisely, the number of words in {acau_red_suff(s)} that are strictly
        smaller than {w} (or strictly greater, if {reverse == TRUE}).
        The result is defined even if {w} is not in {acau_red_suff(s)}.  Cost:
        $O(n(d + 1))$, where $n$ is the length of {w}, and $d$ is the
        average value of {acau_out_deg(t)} over all states {t} in the path
        spelled by {w}. */

      acau_red_add_string(s: acau_red_state_t; READONLY w: String): acau_red_state_t RAISES {acau_red_error_FULL};
      /*
        Returns the state {s'} such that {acau_red_suff(s') == acau_red_suff(s) + {w}}.
        (In particular, returns {s} itself if {w} is already in {acau_red_suff(s)}.

        Note that {acau_red_add_string} does NOT modify the current root
        state, and does not change {acau_red_pref(t)} for any state {t}. */

      acau_red_sub-string(s: acau_red_state_t; READONLY w: String): acau_red_state_t RAISES {acau_red_error_FULL};
      /*
        Returns the state {s'} such that {acau_red_suff(s') == acau_red_suff(s) - {w}}.
        In particular, returns {s} itself if {w} is not in acau_red_suff(s).

        Note that acau_red_sub-string does NOT modify the current root
        state, and does not change {acau_red_pref(t)} for any state {t}.

        In the current implementation, the spacetime cost of
        {acau_red_add_string} and {acau_red_sub-string} is $\Theta(m n)$ in the worst case,
        where $m$ is the number of symbols in the automaton's alphabet,
        and $n$ is the length of {w}. */
        
      acau_red_add_string_maybe_crunch(
          acau_red_state_t *s; 
          String READONLY w; 
          keep: REF States = NULL;
        ): acau_red_state_t;
      
      acau_red_sub_string_maybe_crunch(
          acau_red_state_t *s;
          String READONLY w; 
          keep: REF States = NULL;
        ): acau_red_state_t;
      /*
        Like {acau_red_add_string} and {acau_red_sub-string}, except that if there isn't enough
        space in the automaton, calls {acau_red_crunch(s & keep)}, and tries again.
        If {acau_red_crunch} doesn't return enough space, expands the automaton
        and tries again, until it succeeds (or runs out of memory). */
        
      /***************************************************************************/
      /* ENUMERATION                                                             */
      /***************************************************************************/

      acau_red_enum_out_arcs(s: acau_red_state_t; action: acau_arc_action_t) RAISES {acau_disp_SKIP, acau_disp_STOP};
      /*
        Enumerates the arcs out of state {s}, in the intrinsic order
        (from {acau_first(s)} to {acau_last(s)}), and calls {action(i,a)} for each;
        where {i} is the index of the arc (counting from 0) and {a} the
        Arc data. 

        If {action(i, a)} raises {acau_disp_SKIP}, {acau_red_enum_out_arcs} ignores all
        arcs out of {s} that follow {a} in the intrinsic ordering.
        {acau_red_enum_out_arcs} passes to the client any {acau_disp_STOP}s raised by
        {action}.*/

      acau_red_enum_in_arcs(s: acau_red_state_t; action: acau_arc_action_t) RAISES {acau_disp_SKIP, acau_disp_STOP};
      /*
        Enumerates the arcs that enter {s}, and calls {action(i, a)} for
        each arc, where {i} is the index of the arc in the enumeration
        (counting from 0), and {a} is the data for the arc with its
        direction reversed (i.e., {a.dest} is actually the origin of 
        the arc; the actual destination is always {s}).

        The arcs are visited in order of increasing origin state
        {a.dest}; arcs with same origin {t} are visited in intrinsic 
        order, from {acau_first(t)} to {acau_last(t)}.

        The exceptions {acau_disp_SKIP} and {acau_disp_STOP} are handled just as in {acau_red_enum_out_arcs}. */

      acau_enum_paths(
          acau_red_state_t s;
          enter: acau_red_state_action_t = NULL;
          push: acau_red_path_action_t = NULL;
          pop: acau_red_path_action_t = NULL;
          exit: acau_red_state_action_t = NULL;
        ) RAISES {acau_disp_STOP};
      /*
        Enumerates paths starting from {s}, in standard double order.
        Only paths that consist entirely of proper arcs and states are
        considered.
        
        The method actually enumerates only a subset of all such
        paths, as determined by the client-provided procedures {enter,
        }push", {pop}, and {exit}.

        The enumeration algorithm can be described in terms of a conceptual
        {current path} that initially contains just the state {s},
        and grows or shrinks one arc at a time.  {acau_enum_paths} calls

            {enter(len,o,f)} when the path has {len} arcs and has just landed
                   on a state {o}, whose {final} bit is {f}.

            {push(len,o,i,a)} whenever the current path has {len} arcs,
                   ends at state {o}, and is being extended with the
                   arc {a}, which is the {i}th arc out of {o} (where
                   {i == 0} means arc {acau_first(o)}).

            {pop(len,o,i,a)} whenever the current path has {len+1} arcs
                   and acau_enum_paths is about to remove its last arc {a},
                   which is the {i}th arc out of {o}.

            {exit(len,o,f)} whenever the current path has {len} arcs and
                   ends in state {o} with final bit {f}, and no extensions
                   of it remain to be enumerated.

        Any of the four client actions can stop the enumeration by raising {acau_disp_STOP}.

        Unless the enumeration is aborted, every call to {enter} will be followed
        eventually by a matching call to {exit}, and every call to {push}
        will be followed eventually by a matching call to {pop}.  The typical
        action pattern for a generic node {t} with {n} outgoing proper arcs is

        |     enter(len,t,f)
        |       push(len,t,0,a)
        |         ...
        |       pop(len,t,0,a)
        |       push(len,t,1,b)
        |         ...
        |       pop(len,t,1,b)
        |       ...
        |       push(len,t,n-1,x)
        |         ...
        |       pop(len,t,n-1,x)
        |     exit(len,t,f)

        In particular, if the starting state {s} is non-null, {acau_enum_paths} will call
        {enter(0,s,acau_red_is_final(s))} at the very beginning, and (if not aborted)
        will call {exit(0,s,acau_red_is_final(s))} at the very end.
        (If the starting state {s} is {acau_red_state_NULL}, {acau_enum_paths} will do nothing.)

        Note that the same state {o} may be {enter}ed and {exit}ed many, many times.

        The client can prune branches of the path tree by raising the
        {acau_disp_SKIP} exception in the action procedures:

            if {enter(len,o,f)} raises {acau_disp_SKIP}, {acau_enum_paths} will omit the 
            entire branch of the path tree rooted at current path, and call 
            {exit(len,o,f)} right away, as if {o} had no outgoing arcs;

            if {push(len,o,i,a)} raises {acau_disp_SKIP}, {acau_enum_paths} will
            call {pop(len,o,i,a)} right away, as if the arc {a} led
            to {acau_red_state_NULL};

            if {pop(len,o,i,a)} raises {acau_disp_SKIP}, {acau_enum_paths} ignores any remaining
            arcs out of {o}, and calls {exit(len,o,f)} right away.

        (It is a checked error for {exit(len,o,f)} to raise {acau_disp_SKIP}.)

        If any action is not specified, {acau_enum_paths} will provide a
        a trivial default action that does nothing (and raises no exceptions).
        */

      acau_red_enum_strings(s: acau_red_state_t; enter, exit: acau_red_string-action_t = NULL) RAISES {acau_disp_STOP};
      /*
        Enumerates strings spelled by paths out of {s}, in standard double
        order.

        {acau_red_enum_strings} is similar {acau_enum_paths}, but gives to the {enter} and {exit}
        actions the whole string spelled by the current path, and its last
        state.

        Specifically, {acau_red_enum_strings} calls

            {enter(w, d, f)} when it has reached for the first time
                some state {d}, whose {final} bit is {f}, by spelling string {w},
                and it is about to enumerate all suffixes of {s} whose 
                proper prefix is {w};

            {exit(w, d, f)} when it is going to leave state {d}, whose
                {final} bit {f}, and has just finished enumerating
                all suffixes of {s} whose proper prefix is {w}.

        Either {enter} or {exit} can stop the enumeration by raising {acau_disp_STOP}.

        Unless the enumeration is aborted, every call to {enter} will be followed
        eventually by a matching call to {exit}.  The very first call to {enter}
        will be {enter({}, s, f)}, and the very last call to {exit} will
        be {exit({}, s, f)}, where {{}} is the empty {String}.

        If {enter(w, d, f)} raises {acau_disp_SKIP}, {acau_red_enum_strings} will omit the entire branch
        of the path tree rooted at the resulting current path, and will call
        {exit(w, d, f)} right away.   (It is an error for {exit} to raise {acau_disp_SKIP}.)

        Note that {acau_red_enum_strings} enumerates all the strings that lead from {s}
        to *any* proper state, not just those that lead to a final state.
        If you only want the latter, use {acau_red_enum_suffs}.

        If either {enter} or {exit} is not specified,  {acau_red_enum_strings} supplies
        a trivial default {acau_red_string-action_t} that does nothing (and raises no exceptions).
        */

      acau_red_enum_states(
          States READONLY base;
          substates: BOOL = FALSE;
          enter: acau_red_state_action_t = NULL;
          exit: acau_red_state_action_t = NULL;
        ) RAISES {acau_disp_STOP};
      /*
        The procedure enumerates all states reachable from the states
        base[i], and calls the client procedures {enter} and {exit}
        exactly once each for every non-null state it sees.

        If {substates == FALSE} (the default), the procedure only
        enumerates those states that can be reached from {base} by
        zero or more {R} steps. If {substates == TRUE}, the procedure
        also enumerates all substates of those states; in other words,
        all states that can be reached from {base} by
        zero or more {R} or {acau_rest} operations.
        In either case, {acau_red_enum_states} visits each state exactly once,
        even if it is a substate shared by two or more reachable states.

        The states are enumerated twice in numerical order (which is consistent
        with both the {R} and {acau_rest} operations): once in decreasing
        order with {enter}, and once in increasing order with {exit}.

        More precisely, {acau_red_enum_states} will call {enter(len,s,f)} for each reachable
        state {s}, in *decreasing* numerical order, where {len} is the minimum
        number of {R} steps (or {R} and {acau_rest} steps, if {substates==TRUE})
        from {base} to {s}, and {f} is {acau_red_is_final(s)}.  Note that this call will happen
        before any successors or substates of {s} are {enter}ed.

        That done, {acau_red_enum_states} will call {exit(len,s,f)} for each reachable
        state {s}, in *increasing* numerical order, where {len} and {f} have the
        same meaning as before.  Note that this call will happen after
        all successors and reachable substates of {s} have been {exit}ed.

        The actions {enter} and {exit} can stop the enumeration
        at any time by raising {acau_disp_STOP}.

        The {enter} procedure can also control the set of visited vertices through
        its returned value. If {enter(len,s,f)} raises {acau_disp_SKIP}, {acau_red_enum_states} will
        ignore any arcs out of {s}.  In that case, the states that can be reached
        from {s} will be visited only if they are reachable through other paths
        that do not include {s}.  (It is an error for {exit} to raise {acau_disp_SKIP}.)

        If either {enter} or {exit} is not specified, {acau_red_enum_states} will
        supply a trivial {acau_red_state_action_t} that does nothing (and raise no exceptions)
   .
        */

      /***************************************************************************/
      /* PREFIXES)  AND  AND  (SUFFIXES                                                   */
      /***************************************************************************/

      acau_red_num_prefs(s: acau_red_state_t): uint32_t;
      /*
        The number of distinct paths that lead from the current
        root state to the state {s}. The result is zero
        if the current root is {acau_red_state_NULL}, if {s} is acau_red_state_NULL,
        or {s} is not reachable from the current root. */

      acau_red_num_suffs(s: acau_red_state_t): uint32_t;
      /*
        The number of distinct paths leading from state {s} to some
        final state. The result is zero iff {s} is {acau_red_state_NULL}. */

      acau_red_num_pref_letters(s: acau_red_state_t): uint32_t;
      /*
        The total length of all prefixes of {s}, that is,
        of all paths leading from the current root state to state {s}.
        The result is zero if the current root is {acau_red_state_NULL}, if {s} 
        is {acau_red_state_NULL}, or {s} is not reachable from the current root. */

      acau_red_num_suff_letters(s: acau_red_state_t): uint32_t;
      /*
        The total length of all suffixes of {s}, that is,
        of all paths leading from {s} to some final state.
        The result is zero iff {s} == acau_red_state_NULL or {s == acau_red_state_UNIT}. */

      acau_red_first_prefix(s: acau_red_state_t; e: acau_encoding.h.T): char *;
      /*
        The first string (in the funny order defined by {acau_red_enum_prefs})
        leading from the root to {s}, converted to a char *by applying
        {lpr} to each symbol.   An error if {s} is {acau_red_state_NULL}, or
        is not reachable from the current root (in particular,
        if the current root is {acau_red_state_NULL}). */

      acau_red_first_suffix(s: acau_red_state_t; e: acau_encoding.h.T): char *;
      /*
        The first string (in standard pre-order) leading from {s}
        to some final state, converted to char *by applying {lpr} to
        each symbol.  An error if {s} is {acau_red_state_NULL}. */

      acau_red_full_label(s: acau_red_state_t; e: acau_encoding.h.T; sep: char *= {:}): char *;
      /*
        Computes a label for the state {s}, of the form
        {<first-prefix><sep><first-suffix>}.
        An error if {s} is not reachable from the root,
        or if {s} is {acau_red_state_NULL}. */

      acau_red_enum_prefs(s: acau_red_state_t; action: acau_red_prefix_action_t) RAISES {acau_disp_STOP};
      /*
        Enumerates {acau_red_pref(s)}, the set of all strings leading from the
        root to state {s}, calling {action} on each.
        The state {s} must not be {acau_red_state_NULL}.

        IMPORTANT (1): each string passed to {action} is actually
        the symbol-reversal of the prefix in question.

        IMPORTANT (2): The strings are enumerated in a funny order
        (depth-first order based on {acau_red_enum_in_arcs}) that is NOT
        lexicographic order, and is not related to the standard
        path orders.

        {acau_red_enum_prefs} will return without calling {action}
        if {s} is not reachable from the current root.

        The {action} may raise {acau_disp_STOP} to stop the enumeration. */

      acau_red_enum_suffs(s: acau_red_state_t; action: acau_red_suffix_action_t) RAISES {acau_disp_STOP};
      /*
        Enumerates {acau_red_suff(s)}, the set of all strings leading from
        {s} to some final state, calling {action} on each string.

        The strings are enumerated in standard path pre-order.

        {acau_red_enum_suffs} will return without calling {action}
        iff {s} is {acau_red_state_NULL}.

        The {action} may raise {acau_disp_STOP} to stop the enumeration. */

      acau_red_print_prefs(s: acau_red_state_t; spr: acau_string_printer.h.T);
      /*
        Prints all prefixes of state {s}, in the funny order used by {acau_red_enum_prefs}.
        Uses the given {acau_string_printer.h.T} for formatting control.
        An error if {s} is {acau_red_state_NULL} or is not reachable from the root. */

      acau_red_print_suffs(s: acau_red_state_t; spr: acau_string_printer.h.T);
      /*
        Prints the suffixes of state {s}, in standard path pre-order.
        Uses the given {acau_string_printer.h.T} for formatting control.
        Prints nothing iff {s} is {acau_red_state_NULL}. */

      /***************************************************************************/
      /* CONSTRUCTION                                                            */
      /***************************************************************************/

      acau_red_build(
          acau_red_next_string_proc_t next;        /* Client input procedure */
          wr: FILE * = NULL;             /* Writer for messages, progress report, etc. */
          e: acau_encoding.h.T = NULL;        /* acau_symbol_t/CHAR encoding for messages */
          reportInterval: POS = 1000; /* Status reporting interval */
          flagRedundant: BOOL = TRUE; /* TRUE to print warnings on redundant ops */
        ) RAISES {acau_disp_STOP};
      /*
        Modifies the automaton by adding and/or deleting given set of strings.
        
        {acau_red_build} calls {next(s, len, add)} repeatedly to get each string.  
        The {next} procedure should return the string's length in {len}, and
        the string itself in {s^[0..len-1]} (expanding s^ if necessary).  
        Also, {next} should set {add = TRUE} for strings to be
        added, and {add = FALSE} for strings to be deleted.  Finally, {next}
        should raise {Done} when there is no next string.

        {acau_red_build} will automatically crunch and/or expand the automaton as necessary.
        
        If {wr} is not NULL, {acau_red_build} will periodically print to {wr} a one-line
        progress report, giving the string and symbol count and the size of
        the automaton.  The report will be printed after processing
        {reportInterval} additional strings, after each {acau_red_crunch} and {acau_red_expand},
        and before returning. 
        
        The acau_encoding.h {e} is used only to print out strings in error messages.
        It defaults to a {Plainacau_encoding.h.New()}. */

      /***************************************************************************/
      /* STATISTICS                                                              */
      /***************************************************************************/

      acau_red_num_states(READONLY base: States): uint32_t;
      /*
        The number of distinct states reachanble from {base} by {R} paths. 
        Cost: O("acau_red_max_state()"). */

      acau_red_num_arcs(READONLY base: States): uint32_t;
      /*
        The number of distinct arcs reachable from the {base} states by 
        {R} paths.  Cost: O("acau_red_max_state()"). */

      acau_red_count(READONLY base: States): acau_red_counts_t;
      /*
        Computes:

          o the number of distinct non-null states, final states, and arcs 
            in the part of {self} that is reachable from the given {base} 
            states by {R} paths.

          o the number of distinct non-null substates reachable from {base}
            by {R}+"acau_rest" paths;

          o the total count and length of all strings accepted by 
            the states in {base}.

        Every reachable state, arc, final, or sub-state is counted
        only once, irrespective of how many {base} elements lead to it.
        On the other hand, if two {base} states accept the same string
        (in particular, if they are the same state), that string
        contributes twice to the string and length counts.
        */

      /***************************************************************************/
      /* STORAGE MANAGEMENT                                                      */
      /***************************************************************************/

      acau_red_max_alloc_state(): POS;
      /*
        Implementation data: maximum number of states/substates
        for which storage has been allocated.

        In the current implementation, each distinct non-null state
        or substate  (reachable or not) occupies one unit of storage.
        Attempts to create more than {acau_red_max_alloc_state()} distinct substates
        will generally raise the {acau_red_error_FULL} exception. */

      acau_red_expand(newSize: uint32_t = 0) RAISES {acau_red_error_FULL};
      /*
        Reallocates enough storage for up to {newSize} new substates.
        Preserves all current states and their numbers.
        If {newSize} is not given, uses a constant factor
        (between 1 and 2) times the current size. Raises {acau_red_error_FULL}
        iff the automaton already contains more than {newSize}
        distinct non-null substates. */

      acau_red_crunch(keep: REF States = NULL);
      /*
        Discards all states not reachable from the current root
        or from the states listed in {keep} (if any); then
        squeezes the reachable ones together as closely as possible,
        and updates the current root and the {keep} vector to reflect
        the new acau_red_state_t numbering.

        WARNING: By definition, {acau_red_crunch} generally changes the numbers
        of all states, and deletes some of them.  The client should
        copy all its important acau_red_state_t variables into the {keep} vector,
        and copy them back after the {acau_red_crunch} is done. */

  ;};

typedef
  acau_arc_action_t == PROCEDURE(i: uint32_t; a: Arc) RAISES {acau_disp_SKIP, acau_disp_STOP};
  acau_red_path_action_t == PROCEDURE(len: uint32_t; org: acau_red_state_t; i: uint32_t; a: Arc) RAISES {acau_disp_SKIP, acau_disp_STOP};
  acau_red_state_action_t == PROCEDURE(len: uint32_t; s: acau_red_state_t; final: BOOL) RAISES {acau_disp_SKIP, acau_disp_STOP};
  acau_red_string-action_t == PROCEDURE(READONLY w: String; s: acau_red_state_t; final: BOOL) RAISES {acau_disp_SKIP,acau_disp_STOP};
  acau_red_prefix_action_t == PROCEDURE(READONLY w: String) RAISES {acau_disp_STOP};
  acau_red_suffix_action_t == PROCEDURE(READONLY w: String) RAISES {acau_disp_STOP};

  acau_red_next_string_proc_t == PROCEDURE(
    VAR /*IO*/ s: REF String;  /* The string buffer */
    VAR /*OUT*/ len: uint32_t;      /* Length of returned string */
    VAR /*OUT*/ add: BOOLEAN;  /* TRUE to add, FALSE to delete */
  ) RAISES {Done, acau_disp_STOP};

typedef
  acau_red_counts_t == RECORD
      uint32_t strings;   /* Total strings accepted from all roots */
      uint32_t symbols;   /* Total length os all strings accepted from all roots */
      uint32_t states;    /* Distinct non-null states reachable by {R} */
      uint32_t finals;    /* Distinct reachable states that are final */
      uint32_t arcs;      /* Total number of outgoing arcs from all reachable states */
      uint32_t substates; /* Distinct non-null substates reachable by {R}+"acau_rest" */
    ;};

/***************************************************************************/
/* CREATION                                                                */
/***************************************************************************/

PROCEDURE acau_red_new(size: POS = 1): T;
/*
  Creates a new reduced automaton, initially with only the
  {acau_red_state_NULL} and {acau_red_state_UNIT}, and whose root is {acau_red_state_NULL}.

  The automaton will have space for {size} distinct non-null substates. */

PROCEDURE acau_red_copy(
    T aut;
    size: POS = 1;
    VAR /*OUT*/ map: REF States;
  ): T RAISES {acau_red_error_FULL};
/*
  Returns a new automaton that contains {acau_red_state_NULL}, {acau_red_state_UNIT}, and an isomorphic
  copy of the part of {aut} that is reachable from {aut.acau_red_root()} by chains of
  {R} or {acau_rest} (i.e. all reachable states of {aut}, and their substates).

  The procedure also returns a vector {map} that gives the correspondence
  between old and new state numbers: specifically, {map[s]} is the state of the
  copy that corresponds to state {s} of {aut}.  If {s} has no corresponding
  state in the copy, the preocedure sets {map[s]} to acau_red_state_NULL.

  The {size} parameter has the same meaning as in {New}.
  The procedure may raise {acau_red_error_FULL} if {size} is not enough.
  */

PROCEDURE acau_red_copy_states(
    T from;
    T to;
    acau_red_state_t s;
    VAR /*IO*/ map: REF States;
  ): acau_red_state_t RAISES {acau_red_error_FULL};
/*
  Copies into the {to} automaton all states of the {from} automaton
  that are reachable from {root} by {R} chains.

  If the {map} argument is not NULL , {Copy} assumes that any state {t} with
  {map[t]!=acau_red_state_NULL} has already been copied into the {to} automaton, and its
  number there is {map[t]}.  In that case, {Copy} will set {map[t]}
  appropriately for every state that it copies.  (Note that if {map} is
  non-NULL, one must have {LAST(map^) >= s}. */

/***************************************************************************/
/* IO                                                                      */
/***************************************************************************/

PROCEDURE acau_red_dump(wr: FILE *; aut: T);
PROCEDURE ReducedLoad(rd: Rd.T; minSize: POS = 1): T;
/*
  {acau_dag_dump} writes {aut} to the given writer, in a format that can be
  read back with {Load}.  
  
  The output is in ASCII, but it is not meant to be human-readable. (except for
  the {doc} text, which is written at the beginning of the file.)

  {acau_dag_dump} only writes the existing states and arcs, from {1} to
  {aut.acau_red_max_state()}. The {acau_reduced_t} returned by {Load} is just large enough to
  contain those arcs, or {minSize} arcs, whichever is larger. */

PROCEDURE acau_red_print(wr: FILE *; aut: T; e: acau_encoding.h.T);
/*
  Writes to {wr} a legible description of the automaton states that
  are reachable from the current root state, in the format

    | root == nnnn
    | 
    | 418 *
    |   a -> 4615
    |   b -> 4634
    |   h -> 1950
    | 
    | 4615
    |   a -> 5019
    |   ...

  Each paragraph describes one state. The first line of each paragraph gives
  the state number, and a '*' iff the state is final.  The remaining lines give
  the symbol labels and the destinations of the arcs out of the state.

  The {acau_symbol_t} label of each arc is printed using the procedure
  {e.PrintLetter}. */

PROCEDURE acau_red_print_counts(wr: FILE *; READONLY ct: acau_red_counts_t);
/*
  Prints the given counts to {wr}. */

;} Reduced.
