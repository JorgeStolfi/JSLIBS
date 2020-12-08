INTERFACE Reduced;

(* Representation of reduced deterministic acyclic automata *)
(* See the copyright and disclaimer note at the end of this file. *)

(*
  A "Reduced.T" is a deterministic acyclic finite automaton,
  reduced but not necessarily connected.

  STATES, ARCS: a "Reduced.T" is a directed graph with a finite number
  of vertices, the ``proper states''.  Each state has zero or more
  outgoing directed edges, the ``proper arcs'' of that state,
  each labeled with a distinct symbol and leading to another proper state.

  FINAL BIT: each state has a "final" attribute that may be TRUE or FALSE.

  ACYCLICITY: the directed graph defined by the proper states and proper
  arcs of a "Reduced.T" is acyclic, and from every state there is at least
  one path leading to a final state.

  NULL STATE: it is convenient to assume that every "Reduced.T" has,
  in addition to the proper states described above, a "virtual" state
  "NullState" with "final=FALSE" and no arcs.  It is also convenient
  to assume that whenever a state has no proper transition labeled
  with a given symbol, it has a ``virtual'' transition labeled with
  that symbol and leading to "NullState".  In particular, the
  "NullState" has a virtual transition to itself under any symbol.
  (Thus, the graph is NOT acyclic if "NullState" and it virtual arcs
  are included.)

  ROOT STATE: A "Reduced.T" has a distinguished ``root state'' "Root()", which
  affects the meaning of certain operations (notably those having to do
  with prefix strings, see below).

  In general, there may be states in the automaton that are not reachable
  from its current root.  In fact, there may be no single state from
  which all other states are reachable.  In particular, the root state
  may be "NullState".

  SUFFIXES OF A STATE: each proper state "s" is said to "accept" or "recognize"
  a language "Suff(s)", the set of ``suffixes'' of "s", which is the set of all
  "String"s spelled by directed paths that lead from "s" to some final state.
  Obviously, the the language "Suff(s)" for a proper state is always finite
  and non-empty.  By definition, "Suff(NullState)" is the empty set.

  PREFIXES OF A STATE: each proper state "s" has also a set of ``prefixes''
  "Pref(s)", the set of all "String"s spelled by paths that lead from the
  root state to "s".  In particular, if "s" is the current root state,
  then "Pref(s)" contains only the empty string; and if "s" is not reachable
  from the root state, then "Pref(s)" is the empty set.
  
  For technical and logical reasons, the set "Pref(NullState)" is undefined.

  STATE UNIQUENESS: the implementation ensures that the automaton
  is ``reduced'', that is, no two distinct states (proper or null)
  have the same set of suffix strings.  In particular, the "NullState" is
  the only state that accepts the empty set of strings {}.

  DETERMINISM: according to the definitions above, the outgoung arcs
  of any given state are labeled with distinct symbols.  It follows
  that the automaton is deterministic: given any state "s" and any
  String "w", there is exactly one state "t" (possibly null) that
  can be reached from "s" by a path whole labels spell "w".
  Therefore, any two distinct states have disjoint prefix sets;
  and any proper reachable state can be uniquely identified by any
  of its prefix strings.

  UNIT STATE: besides the "NullState", every "Reduced.T" has
  exactly one proper state, the "UnitState", that has "final=TRUE" and no
  outgoing arcs.  Obviously, that state accepts the language consisting
  only of the empty string {()}.
  
  INTRINSIC ARC ORDER: for every state "s" that has outgoing arcs
  (that is, any state other than "NullState" or "UnitState"), there is
  one distinguished arc "Last(s)" out of "s", and another state
  "Rest(s)", such that the set of arcs leaving the latter are those
  leaving "s", minus the arc "Last(s)".
  
  The methods "Last" and "Rest" define the ``intrinsic ordering'' for the
  arcs out "s".  Note that "Rest" moves {\em backwards} in this ordering.
    
  SUBSTATES: A "substate" of a state "s" of a Reduced.T is a state "t"
  that whose transitions are an initial prefix of the transitions out
  of "s", in intrinsic order (or, in other words, that can be obtained
  from "s" by zero or more applications of the "Rest" method).

  In particular, the "NullState" is a substate of every state
  with "final=FALSE", and the "UnitState" is a substate of every
  state with "final=TRUE".

  If a "Reduced.T" includes a state "s", it automatically includes
  all substates of "s".
  
  TOPOLOGICAL ORDERING OF THE STATES: The states in a "Reduced.T" are
  identified by numbers, assigned as they are created.  The numerical
  order of the states is compatible with the transitions and the
  substate relation; that is, every state is greater (numerically)
  than any of its proper substates, and also greater than the
  destination of any of its outgoing arcs.

  STANDARD PATH ORDERS: The intrinsic ordering of the arcs out of each
  state induces a ``standard'' order for the set of all paths out of a
  fixed state "s".  In this order, all paths that leave "s" through
  the arc "Last(s)" occur after those that leave through other arcs.
  
  There are actually three standard path orderings, differing on where
  the paths that {\em end} at "s" appear: either before all the paths
  that {\em go through} "s" (``pre-order''), or after them
  (``end-order''), or both before and after them (``double order'').
  
  STRING SPELLING: Every path out of a state "s" is said to
  ``spell out'' the "String" consisting of the symbols that label the
  arcs along the path.  Conversely, since the automaton is
  deterministic, each "String" defines a unique path out of a given
  state "s".
  
  LEXICOGRAPHIC ORDER: The ordering of the "Symbol" values induces a
  ``lexicographic ordering'' of any finite set of "String"s, and hence
  of all the paths out of a state "s"
  
  Note that the lexicographic ordering is distinct from all the
  standard path orderings (unless the arcs in every state "s" happen
  to be sorted by symbol value, with "Last(s)" being the highest).

  IMMUTABILITY:  Except for the "Crunch" method (below), the procedures
  and methods in this interface never modify or destroy existing states,
  but merely create new ones as needed.  Therefore, if "Crunch" is not used,
  the "final" bit and outgoing arcs of a state "s" are constant
  throughout its lifetime, and so are any properties of "s"
  that can defined in terms of these (such as the suffix set).

  Moreover, "Crunch", "SetRoot", and "Build" are the only two methods or
  procedures in this interface that change the current root state.  It follows
  the prefix set Pref(s) of any state "s" is constant as long as these methods
  are not used.  *)

IMPORT Basics, Rd, Wr, StringPrinter, Encoding;
FROM Basics IMPORT NAT, POS, BOOL;
FROM Basics IMPORT Done, Full, Skip, Abort;

(*
  The methods and procedures below raise "Basics.Full" if there isn't
  enough allocated storage in the automaton for the requested operation. *)

TYPE
  Symbol = [ FIRST(Basics.Symbol)+1 .. LAST(Basics.Symbol) ];
  (*
    Note that FIRST(Basics.Symbol) is reserved for internal
    use of the implementation *)

  String = Basics.String;
  (*
    Should be ARRAY OF Reduced.Symbol, excluding FIRST(Basics.Symbol),
    but that would require special versions of ExpandString, etc. *)

  State = NAT;
  States = ARRAY OF State;
  
  Arc = RECORD symbol: Symbol; dest: State END;
  Arcs = ARRAY OF Arc;

CONST
  NullState = 0; (* The unique (improper) state with no arcs and "final=FALSE" *)
  UnitState = 1; (* The unique (proper) state with no arcs and "final=TRUE" *)

TYPE
  T <: Public;

  Public = OBJECT
  
      doc: TEXT;  
        (* 
          An arbitrary multi-line text that describes the automaton. 
          If not empty, it should end with a newline character. *)
      
    METHODS

      (*
        Note: the methods that deal with prefixes ("NPrefs", "EnumPrefs", ...)  and
        those that deal with incoming edges ("InDeg", "EnumInArcs") automatically
        construct internal tables whose cost is proportional to "MaxState()".
        The tables are reused by subsequent calls to those methods, but are
        invalidated by "Crunch" and "SetRoot".  So, the first prefix-related
        method call after a "Crunch" or "SetRoot" may be a lot more expensive
        than ordinary calls. *)

      MaxState(): POS;
      (*
        The number of the highest-numbered state in the automaton.
        Note: not all numbers in "[0..MaxState()]" correspond
        to currently reachable states.
        Cost: $O(1)$ time, 0 space. *)

      Final(s: State): BOOL;
      (*
        TRUE iff "s" is a final state of the automaton.
        Cost: $O("OutDeg(s)")$ time, 0 space. *)

      HasArcs(s: State): BOOL;
      (*
        TRUE iff the state "s" has at least one proper outgoing arc.
        Cost: $O(1)$ time, 0 space. *)

      Last(s: State): Arc;
      (*
        The last proper arc out of "s", in intrinsic order.
        Requires "HasArcs(s)".
        Cost: $O(1)$ time, 0 space. *)

      Rest(s: State): State;
      (*
        The state obtained from "s" by removing the transition "Last(s)",
        with the same final bit as "s". Requires "HasArcs(s)".
        Cost: $O(1)$ time, 0 space. *)

      Root(): State;
      (*
        Returns the current root state.
        Cost: $O(1)$ time, 0 space. *)

      SetRoot(s: State);
      (*
        Makes "s" the current root state.  Cost: $O(1)$ time, 0 space.
        However, note that this operation may implicitly change the
        prefix sets of many states. *)

      (******************************************************************)
      (* DERIVED METHODS                                                *)
      (******************************************************************)

      (*
        These methods could be ordinary procedures, since
        they can be defined in terms of the primitives above;
        they are declared here as methods to reduce client confusion. (?) *)

      OutDeg(s: State): NAT;
      (*
        Number of proper arcs out of state "s" (0 iff "s = NullState"
        or "s = UnitState"). Cost: $O("OutDeg(s)")$ time, 0 space. *)

      InDeg(s: State): NAT;
      (*
        Number of proper arcs that lead to state "s" from some state reachable
        from "Root()".  Returns 0 is "s = Root()" or "s = NullState".
        Cost: $O("InDeg(s)")$ time, 0 space (plus the cost of
        rebuilding the reverse-data tables, if "Root()" has changed). *)

      First(s: State): Arc;
      (*
        The first arc out of state "s", in intrinsic order. 
        Requires "s # NullState". Cost: $O("OutDeg(s)")$ time, 0 space. *)

      R(s: State; symbol: Symbol): State;
      (*
        The successor of the given state through the arc labeled with
        the given symbol.  Returns "NullState" if "s" is "NullState" or has
        no outgoing edge labeled with that "symbol".
        Cost: $O("OutDeg(s)")$ time, 0 space. *)

      SetArc(s: State; symbol: Symbol; dest: State): State RAISES {Full};
      (*
        The unique state that has the same arcs and final bit as "s",
        except that it goes to "dest" through the given symbol.
        Creates the state if it doesn't exist yet.

        In particular, returns "s" if "R(s, symbol) = dest"
        already.  The "dest" state may be "NullState", in which case
        the arc from "s" through the given symbol is removed.
        Cost: $O("OutDeg(s)")$ time and space, if hashing works. *)

      SetFinal(s: State; final: BOOLEAN): State RAISES {Full};
      (*
        The unique state that has the same arcs as "s", except that
        its final bit is "final".  Creates the state if it doesn't
        exist yet.

        The resulting state accepts all the strings accepted by "s",
        plus the empty string.  In particular, the call
        "SetFinal(NullState, TRUE)" returns "UnitState".

        Cost: $O("OutDeg(s)")$ time and space, if hashing works.
        Hence, when building a state from scratch, it is better to set
        the final bit before adding the outgoing edges. *)

      Walk(s: State; READONLY w: String): State;
      (*
        Spells "w" in the automaton starting from "s", that is,
        returns the state reached from "s" by the unique path
        whose edges are labeled with the symbols of "w". *)

      Accepts(s: State; READONLY w: String): BOOL;
      (*
        TRUE iff "w" is in the language "Suff(s)".  Equivalent to
        "Final(Walk(s, w))". *)

      Rank(s: State; READONLY w: String; reverse: BOOL := FALSE): NAT;
      (*
        The rank of "w" in "Suff(s)", in lexicographic order.  More
        precisely, the number of words in "Suff(s)" that are strictly
        smaller than "w" (or strictly greater, if "reverse = TRUE").
        The result is defined even if "w" is not in "Suff(s)".  Cost:
        $O(n(d + 1))$, where $n$ is the length of "w", and $d$ is the
        average value of "OutDeg(t)" over all states "t" in the path
        spelled by "w". *)

      AddString(s: State; READONLY w: String): State RAISES {Full};
      (*
        Returns the state "s'" such that "Suff(s') = Suff(s) + {w}".
        (In particular, returns "s" itself if "w" is already in "Suff(s)".

        Note that "AddString" does NOT modify the current root
        state, and does not change "Pref(t)" for any state "t". *)

      SubString(s: State; READONLY w: String): State RAISES {Full};
      (*
        Returns the state "s'" such that "Suff(s') = Suff(s) - {w}".
        In particular, returns "s" itself if "w" is not in Suff(s).

        Note that SubString does NOT modify the current root
        state, and does not change "Pref(t)" for any state "t".

        In the current implementation, the spacetime cost of
        "AddString" and "SubString" is $\Theta(m n)$ in the worst case,
        where $m$ is the number of symbols in the automaton's alphabet,
        and $n$ is the length of "w". *)
        
      AddStringMaybeCrunch(
          VAR s: State; 
          READONLY w: String; 
          keep: REF States := NIL;
        ): State;
      
      SubStringMaybeCrunch(
          VAR s: State;
          READONLY w: String; 
          keep: REF States := NIL;
        ): State;
      (*
        Like "AddString" and "SubString", except that if there isn't enough
        space in the automaton, calls "Crunch(s & keep)", and tries again.
        If "Crunch" doesn't return enough space, expands the automaton
        and tries again, until it succeeds (or runs out of memory). *)
        
      (***************************************************************************)
      (* ENUMERATION                                                             *)
      (***************************************************************************)

      EnumOutArcs(s: State; action: ArcAction) RAISES {Skip, Abort};
      (*
        Enumerates the arcs out of state "s", in the intrinsic order
        (from "First(s)" to "Last(s)"), and calls "action(i,a)" for each;
        where "i" is the index of the arc (counting from 0) and "a" the
        Arc data. 

        If "action(i, a)" raises "Skip", "EnumOutArcs" ignores all
        arcs out of "s" that follow "a" in the intrinsic ordering.
        "EnumOutArcs" passes to the client any "Abort"s raised by
        "action".*)

      EnumInArcs(s: State; action: ArcAction) RAISES {Skip, Abort};
      (*
        Enumerates the arcs that enter "s", and calls "action(i, a)" for
        each arc, where "i" is the index of the arc in the enumeration
        (counting from 0), and "a" is the data for the arc with its
        direction reversed (i.e., "a.dest" is actually the origin of 
        the arc; the actual destination is always "s").

        The arcs are visited in order of increasing origin state
        "a.dest"; arcs with same origin "t" are visited in intrinsic 
        order, from "First(t)" to "Last(t)".

        The exceptions "Skip" and "Abort" are handled just as in "EnumOutArcs". *)

      EnumPaths(
          s: State;
          enter: StateAction := NIL;
          push: PathAction := NIL;
          pop: PathAction := NIL;
          exit: StateAction := NIL;
        ) RAISES {Abort};
      (*
        Enumerates paths starting from "s", in standard double order.
        Only paths that consist entirely of proper arcs and states are
        considered.
        
        The method actually enumerates only a subset of all such
        paths, as determined by the client-provided procedures "enter,
        "push", "pop", and "exit".

        The enumeration algorithm can be described in terms of a conceptual
        "current path" that initially contains just the state "s",
        and grows or shrinks one arc at a time.  "EnumPaths" calls

            "enter(len,o,f)" when the path has "len" arcs and has just landed
                   on a state "o", whose "final" bit is "f".

            "push(len,o,i,a)" whenever the current path has "len" arcs,
                   ends at state "o", and is being extended with the
                   arc "a", which is the "i"th arc out of "o" (where
                   "i = 0" means arc "First(o)").

            "pop(len,o,i,a)" whenever the current path has "len+1" arcs
                   and EnumPaths is about to remove its last arc "a",
                   which is the "i"th arc out of "o".

            "exit(len,o,f)" whenever the current path has "len" arcs and
                   ends in state "o" with final bit "f", and no extensions
                   of it remain to be enumerated.

        Any of the four client actions can stop the enumeration by raising "Abort".

        Unless the enumeration is aborted, every call to "enter" will be followed
        eventually by a matching call to "exit", and every call to "push"
        will be followed eventually by a matching call to "pop".  The typical
        action pattern for a generic node "t" with "n" outgoing proper arcs is

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

        In particular, if the starting state "s" is non-null, "EnumPaths" will call
        "enter(0,s,Final(s))" at the very beginning, and (if not aborted)
        will call "exit(0,s,Final(s))" at the very end.
        (If the starting state "s" is "NullState", "EnumPaths" will do nothing.)

        Note that the same state "o" may be "enter"ed and "exit"ed many, many times.

        The client can prune branches of the path tree by raising the
        "Skip" exception in the action procedures:

            if "enter(len,o,f)" raises "Skip", "EnumPaths" will omit the 
            entire branch of the path tree rooted at current path, and call 
            "exit(len,o,f)" right away, as if "o" had no outgoing arcs;

            if "push(len,o,i,a)" raises "Skip", "EnumPaths" will
            call "pop(len,o,i,a)" right away, as if the arc "a" led
            to "NullState";

            if "pop(len,o,i,a)" raises "Skip", "EnumPaths" ignores any remaining
            arcs out of "o", and calls "exit(len,o,f)" right away.

        (It is a checked error for "exit(len,o,f)" to raise "Skip".)

        If any action is not specified, "EnumPaths" will provide a
        a trivial default action that does nothing (and raises no exceptions).
        *)

      EnumStrings(s: State; enter, exit: StringAction := NIL) RAISES {Abort};
      (*
        Enumerates strings spelled by paths out of "s", in standard double
        order.

        "EnumStrings" is similar "EnumPaths", but gives to the "enter" and "exit"
        actions the whole string spelled by the current path, and its last
        state.

        Specifically, "EnumStrings" calls

            "enter(w, d, f)" when it has reached for the first time
                some state "d", whose "final" bit is "f", by spelling string "w",
                and it is about to enumerate all suffixes of "s" whose 
                proper prefix is "w";

            "exit(w, d, f)" when it is going to leave state "d", whose
                "final" bit "f", and has just finished enumerating
                all suffixes of "s" whose proper prefix is "w".

        Either "enter" or "exit" can stop the enumeration by raising "Abort".

        Unless the enumeration is aborted, every call to "enter" will be followed
        eventually by a matching call to "exit".  The very first call to "enter"
        will be "enter({}, s, f)", and the very last call to "exit" will
        be "exit({}, s, f)", where "{}" is the empty "String".

        If "enter(w, d, f)" raises "Skip", "EnumStrings" will omit the entire branch
        of the path tree rooted at the resulting current path, and will call
        "exit(w, d, f)" right away.   (It is an error for "exit" to raise "Skip".)

        Note that "EnumStrings" enumerates all the strings that lead from "s"
        to *any* proper state, not just those that lead to a final state.
        If you only want the latter, use "EnumSuffs".

        If either "enter" or "exit" is not specified,  "EnumStrings" supplies
        a trivial default "StringAction" that does nothing (and raises no exceptions).
        *)

      EnumStates(
          READONLY base: States;
          substates: BOOL := FALSE;
          enter: StateAction := NIL;
          exit: StateAction := NIL;
        ) RAISES {Abort};
      (*
        The procedure enumerates all states reachable from the states
        base[i], and calls the client procedures "enter" and "exit"
        exactly once each for every non-null state it sees.

        If "substates = FALSE" (the default), the procedure only
        enumerates those states that can be reached from "base" by
        zero or more "R" steps. If "substates = TRUE", the procedure
        also enumerates all substates of those states; in other words,
        all states that can be reached from "base" by
        zero or more "R" or "Rest" operations.
        In either case, "EnumStates" visits each state exactly once,
        even if it is a substate shared by two or more reachable states.

        The states are enumerated twice in numerical order (which is consistent
        with both the "R" and "Rest" operations): once in decreasing
        order with "enter", and once in increasing order with "exit".

        More precisely, "EnumStates" will call "enter(len,s,f)" for each reachable
        state "s", in *decreasing* numerical order, where "len" is the minimum
        number of "R" steps (or "R" and "Rest" steps, if "substates=TRUE")
        from "base" to "s", and "f" is "Final(s)".  Note that this call will happen
        before any successors or substates of "s" are "enter"ed.

        That done, "EnumStates" will call "exit(len,s,f)" for each reachable
        state "s", in *increasing* numerical order, where "len" and "f" have the
        same meaning as before.  Note that this call will happen after
        all successors and reachable substates of "s" have been "exit"ed.

        The actions "enter" and "exit" can stop the enumeration
        at any time by raising "Abort".

        The "enter" procedure can also control the set of visited vertices through
        its returned value. If "enter(len,s,f)" raises "Skip", "EnumStates" will
        ignore any arcs out of "s".  In that case, the states that can be reached
        from "s" will be visited only if they are reachable through other paths
        that do not include "s".  (It is an error for "exit" to raise "Skip".)

        If either "enter" or "exit" is not specified, "EnumStates" will
        supply a trivial "StateAction" that does nothing (and raise no exceptions)
   .
        *)

      (***************************************************************************)
      (* PREFIXES AND SUFFIXES                                                   *)
      (***************************************************************************)

      NPrefs(s: State): NAT;
      (*
        The number of distinct paths that lead from the current
        root state to the state "s". The result is zero
        if the current root is "NullState", if "s" is NullState,
        or "s" is not reachable from the current root. *)

      NSuffs(s: State): NAT;
      (*
        The number of distinct paths leading from state "s" to some
        final state. The result is zero iff "s" is "NullState". *)

      NPrefLetters(s: State): NAT;
      (*
        The total length of all prefixes of "s", that is,
        of all paths leading from the current root state to state "s".
        The result is zero if the current root is "NullState", if "s" 
        is "NullState", or "s" is not reachable from the current root. *)

      NSuffLetters(s: State): NAT;
      (*
        The total length of all suffixes of "s", that is,
        of all paths leading from "s" to some final state.
        The result is zero iff "s" = NullState or "s = UnitState". *)

      FirstPrefix(s: State; e: Encoding.T): TEXT;
      (*
        The first string (in the funny order defined by "EnumPrefs")
        leading from the root to "s", converted to a TEXT by applying
        "lpr" to each symbol.   An error if "s" is "NullState", or
        is not reachable from the current root (in particular,
        if the current root is "NullState"). *)

      FirstSuffix(s: State; e: Encoding.T): TEXT;
      (*
        The first string (in standard pre-order) leading from "s"
        to some final state, converted to TEXT by applying "lpr" to
        each symbol.  An error if "s" is "NullState". *)

      FullLabel(s: State; e: Encoding.T; sep: TEXT := ":"): TEXT;
      (*
        Computes a label for the state "s", of the form
        "<first-prefix><sep><first-suffix>".
        An error if "s" is not reachable from the root,
        or if "s" is "NullState". *)

      EnumPrefs(s: State; action: PrefixAction) RAISES {Abort};
      (*
        Enumerates "Pref(s)", the set of all strings leading from the
        root to state "s", calling "action" on each.
        The state "s" must not be "NullState".

        IMPORTANT (1): each string passed to "action" is actually
        the symbol-reversal of the prefix in question.

        IMPORTANT (2): The strings are enumerated in a funny order
        (depth-first order based on "EnumInArcs") that is NOT
        lexicographic order, and is not related to the standard
        path orders.

        "EnumPrefs" will return without calling "action"
        if "s" is not reachable from the current root.

        The "action" may raise "Abort" to stop the enumeration. *)

      EnumSuffs(s: State; action: SuffixAction) RAISES {Abort};
      (*
        Enumerates "Suff(s)", the set of all strings leading from
        "s" to some final state, calling "action" on each string.

        The strings are enumerated in standard path pre-order.

        "EnumSuffs" will return without calling "action"
        iff "s" is "NullState".

        The "action" may raise "Abort" to stop the enumeration. *)

      PrintPrefs(s: State; spr: StringPrinter.T);
      (*
        Prints all prefixes of state "s", in the funny order used by "EnumPrefs".
        Uses the given "StringPrinter.T" for formatting control.
        An error if "s" is "NullState" or is not reachable from the root. *)

      PrintSuffs(s: State; spr: StringPrinter.T);
      (*
        Prints the suffixes of state "s", in standard path pre-order.
        Uses the given "StringPrinter.T" for formatting control.
        Prints nothing iff "s" is "NullState". *)

      (***************************************************************************)
      (* CONSTRUCTION                                                            *)
      (***************************************************************************)

      Build(
          next: NextStringProc;        (* Client input procedure *)
          wr: Wr.T := NIL;             (* Writer for messages, progress report, etc. *)
          e: Encoding.T := NIL;        (* Symbol/CHAR encoding for messages *)
          reportInterval: POS := 1000; (* Status reporting interval *)
          flagRedundant: BOOL := TRUE; (* TRUE to print warnings on redundant ops *)
        ) RAISES {Abort};
      (*
        Modifies the automaton by adding and/or deleting given set of strings.
        
        "Build" calls "next(s, len, add)" repeatedly to get each string.  
        The "next" procedure should return the string's length in "len", and
        the string itself in "s^[0..len-1]" (expanding s^ if necessary).  
        Also, "next" should set "add := TRUE" for strings to be
        added, and "add := FALSE" for strings to be deleted.  Finally, "next"
        should raise "Done" when there is no next string.

        "Build" will automatically crunch and/or expand the automaton as necessary.
        
        If "wr" is not NIL, "Build" will periodically print to "wr" a one-line
        progress report, giving the string and symbol count and the size of
        the automaton.  The report will be printed after processing
        "reportInterval" additional strings, after each "Crunch" and "Expand",
        and before returning. 
        
        The Encoding "e" is used only to print out strings in error messages.
        It defaults to a "PlainEncoding.New()". *)

      (***************************************************************************)
      (* STATISTICS                                                              *)
      (***************************************************************************)

      NStates(READONLY base: States): NAT;
      (*
        The number of distinct states reachanble from "base" by "R" paths. 
        Cost: O("MaxState()"). *)

      NArcs(READONLY base: States): NAT;
      (*
        The number of distinct arcs reachable from the "base" states by 
        "R" paths.  Cost: O("MaxState()"). *)

      Count(READONLY base: States): Counts;
      (*
        Computes:

          o the number of distinct non-null states, final states, and arcs 
            in the part of "self" that is reachable from the given "base" 
            states by "R" paths.

          o the number of distinct non-null substates reachable from "base"
            by "R"+"Rest" paths;

          o the total count and length of all strings accepted by 
            the states in "base".

        Every reachable state, arc, final, or sub-state is counted
        only once, irrespective of how many "base" elements lead to it.
        On the other hand, if two "base" states accept the same string
        (in particular, if they are the same state), that string
        contributes twice to the string and length counts.
        *)

      (***************************************************************************)
      (* STORAGE MANAGEMENT                                                      *)
      (***************************************************************************)

      MaxAllocState(): POS;
      (*
        Implementation data: maximum number of states/substates
        for which storage has been allocated.

        In the current implementation, each distinct non-null state
        or substate  (reachable or not) occupies one unit of storage.
        Attempts to create more than "MaxAllocState()" distinct substates
        will generally raise the "Full" exception. *)

      Expand(newSize: NAT := 0) RAISES {Full};
      (*
        Reallocates enough storage for up to "newSize" new substates.
        Preserves all current states and their numbers.
        If "newSize" is not given, uses a constant factor
        (between 1 and 2) times the current size. Raises "Full"
        iff the automaton already contains more than "newSize"
        distinct non-null substates. *)

      Crunch(keep: REF States := NIL);
      (*
        Discards all states not reachable from the current root
        or from the states listed in "keep" (if any); then
        squeezes the reachable ones together as closely as possible,
        and updates the current root and the "keep" vector to reflect
        the new State numbering.

        WARNING: By definition, "Crunch" generally changes the numbers
        of all states, and deletes some of them.  The client should
        copy all its important State variables into the "keep" vector,
        and copy them back after the "Crunch" is done. *)

  END;

TYPE
  ArcAction = PROCEDURE(i: NAT; a: Arc) RAISES {Skip, Abort};
  PathAction = PROCEDURE(len: NAT; org: State; i: NAT; a: Arc) RAISES {Skip, Abort};
  StateAction = PROCEDURE(len: NAT; s: State; final: BOOL) RAISES {Skip, Abort};
  StringAction = PROCEDURE(READONLY w: String; s: State; final: BOOL) RAISES {Skip,Abort};
  PrefixAction = PROCEDURE(READONLY w: String) RAISES {Abort};
  SuffixAction = PROCEDURE(READONLY w: String) RAISES {Abort};

  NextStringProc = PROCEDURE(
    VAR (*IO*) s: REF String;  (* The string buffer *)
    VAR (*OUT*) len: NAT;      (* Length of returned string *)
    VAR (*OUT*) add: BOOLEAN;  (* TRUE to add, FALSE to delete *)
  ) RAISES {Done, Abort};

TYPE
  Counts = RECORD
      strings: NAT;   (* Total strings accepted from all roots *)
      symbols: NAT;   (* Total length os all strings accepted from all roots *)
      states: NAT;    (* Distinct non-null states reachable by "R" *)
      finals: NAT;    (* Distinct reachable states that are final *)
      arcs: NAT;      (* Total number of outgoing arcs from all reachable states *)
      substates: NAT; (* Distinct non-null substates reachable by "R"+"Rest" *)
    END;

(***************************************************************************)
(* CREATION                                                                *)
(***************************************************************************)

PROCEDURE New(size: POS := 1): T;
(*
  Creates a new reduced automaton, initially with only the
  "NullState" and "UnitState", and whose root is "NullState".

  The automaton will have space for "size" distinct non-null substates. *)

PROCEDURE Copy(
    aut: T;
    size: POS := 1;
    VAR (*OUT*) map: REF States;
  ): T RAISES {Full};
(*
  Returns a new automaton that contains "NullState", "UnitState", and an isomorphic
  copy of the part of "aut" that is reachable from "aut.Root()" by chains of
  "R" or "Rest" (i.e. all reachable states of "aut", and their substates).

  The procedure also returns a vector "map" that gives the correspondence
  between old and new state numbers: specifically, "map[s]" is the state of the
  copy that corresponds to state "s" of "aut".  If "s" has no corresponding
  state in the copy, the preocedure sets "map[s]" to NullState.

  The "size" parameter has the same meaning as in "New".
  The procedure may raise "Full" if "size" is not enough.
  *)

PROCEDURE CopyStates(
    from: T;
    to: T;
    s: State;
    VAR (*IO*) map: REF States;
  ): State RAISES {Full};
(*
  Copies into the "to" automaton all states of the "from" automaton
  that are reachable from "root" by "R" chains.

  If the "map" argument is not NIL , "Copy" assumes that any state "t" with
  "map[t] # NullState" has already been copied into the "to" automaton, and its
  number there is "map[t]".  In that case, "Copy" will set "map[t]"
  appropriately for every state that it copies.  (Note that if "map" is
  non-NIL, one must have "LAST(map^) >= s". *)

(***************************************************************************)
(* IO                                                                      *)
(***************************************************************************)

PROCEDURE Dump(wr: Wr.T; aut: T);
PROCEDURE Load(rd: Rd.T; minSize: POS := 1): T;
(*
  "Dump" writes "aut" to the given writer, in a format that can be
  read back with "Load".  
  
  The output is in ASCII, but it is not meant to be human-readable. (except for
  the "doc" text, which is written at the beginning of the file.)

  "Dump" only writes the existing states and arcs, from "1" to
  "aut.MaxState()". The "Reduced.T" returned by "Load" is just large enough to
  contain those arcs, or "minSize" arcs, whichever is larger. *)

PROCEDURE Print(wr: Wr.T; aut: T; e: Encoding.T);
(*
  Writes to "wr" a legible description of the automaton states that
  are reachable from the current root state, in the format

    | root = nnnn
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

  The "Symbol" label of each arc is printed using the procedure
  "e.PrintLetter". *)

PROCEDURE PrintCounts(wr: Wr.T; READONLY ct: Counts);
(*
  Prints the given counts to "wr". *)

END Reduced.

(****************************************************************************)
(* (C) Copyright 1992 Universidade Estadual de Campinas (UNICAMP)           *)
(*                    Campinas, SP, Brazil                                  *)
(*                                                                          *)
(* Authors:                                                                 *)
(*                                                                          *)
(*   Tomasz Kowaltowski  - CS Dept, UNICAMP <tomasz@dcc.unicamp.br>         *)
(*   Claudio L. Lucchesi - CS Dept, UNICAMP <lucchesi@dcc.unicamp.br>       *)
(*   Jorge Stolfi        - CS Dept, UNICAMP <stolfi@dcc.unicamp.br>         *)
(*                                                                          *)
(* This file can be freely distributed, modified, and used for any          *)
(*   non-commercial purpose, provided that this copyright and authorship    *)
(*   notice be included in any copy or derived version of this file.        *)
(*                                                                          *)
(* DISCLAIMER: This software is offered ``as is'', without any guarantee    *)
(*   as to fitness for any particular purpose.  Neither the copyright       *)
(*   holder nor the authors or their employers can be held responsible for  *)
(*   any damages that may result from its use.                              *)
(****************************************************************************)
