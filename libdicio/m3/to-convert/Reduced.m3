MODULE Reduced;

(* See the copyright and disclaimer note at the end of this file. *)

IMPORT Rd, Wr, TextWr, TextRd, Fmt, Text, Thread;
IMPORT Basics, DAG, StringPrinter, Encoding, PlainEncoding;
FROM Basics IMPORT NAT, POS, BOOL;
FROM Basics IMPORT Done, Full, Skip, Abort;
FROM Basics IMPORT WrNat, ReadParam;

REVEAL
  T = Public BRANDED OBJECT

      dag: DAG.T;
      (*
        The DAG.T that is being used to implement this
        reduced automaton ("aut").

        Each state of "aut" is represented by a state
        of "dag", whose out-arcs all have distinct "rd"-labels
        and zero "wr"-labels, and are sorted by decreasing "rd"-label.
        A state of "aut" is final iff the corresponding state
        of "dag" has an outgoing arc with "rd=NullLetter"; these
        special arcs all lead to the NullState of "dag".

        These conventions, plus the "state uniqueness" invariant
        of a DAG.T, ensure (I hope) that the automaton "aut" is reduced
        in the usual sense. *)

      root: State;
      (*
        The current root state. May be NullState. *)

      (*************************************************************************)
      (* PREFIX/SUFFIX DATA                                                    *)
      (*************************************************************************)

      (*
        The folowing fields are auxiliary tables used by methods that deal
        with prefixes and suffixes.

        The "nSuffs" "nSufLetters" arrays are always allocated, big enough, 
        and up-to-date.

        The prefix-related fields ("rdag", "rev", "dir", "nPrefs", and 
        "nPrefLetters") are (re)computed only when needed (and if "aut.root" 
        is not NullState). If "prefRoot" is NullState, then those fields are
        garbage; if "prefRoot" is a proper state, those fields are valid
        provided that the current root is "prefRoot".

        The reverse DAG "rdag" has one arc from rev[t] to rev[s]
        whenever the automaton has a proper arc from "s" to "t".
        The virtual arcs of "aut" that go to NullState, including
        the invisible arcs labelled with NullLetter that denote final states
        of "aut", have no correspondents in "rdag". On the other hand, "rdag"
        has one arc labelled with "NullLetter" going from rev[prefRoot]
        to NullState. *)

      nSuffs: REF ARRAY OF NAT;       (* Number of suffix strings per state *)
      nSuffLetters: REF ARRAY OF NAT; (* Total length of suffix strings per state *)

      prefRoot: State;     (* Root used to compute "rdag, rev, dir, nPrefs" *)

      rdag: DAG.T;                    (* The DAG with reverse transitions *)
      rev: REF ARRAY OF DAG.State;    (* Maps states of self to states of "rdag" *)
      dir: REF ARRAY OF State;        (* Inverse of "rev", or NullState if undef *)
      nPrefs: REF ARRAY OF NAT;       (* Number of prefix strings for each state *)
      nPrefLetters: REF ARRAY OF NAT; (* Total length of prefix strings per state *)

    OVERRIDES

      MaxState := MaxState;
      Final := Final;
      HasArcs := HasArcs;
      Last := Last;
      Append := Append;
      Rest := Rest;
      Root := Root;
      SetRoot := SetRoot;

      OutDeg := OutDeg;
      InDeg := InDeg;
      First := First;
      R := R;
      SetArc := SetArc;
      SetFinal := SetFinal;
      Walk := Walk;
      Accepts := Accepts;
      Rank := Rank;
      AddString := AddString;
      SubString := SubString;
      AddStringMaybeCrunch := AddStringMaybeCrunch;
      SubStringMaybeCrunch := SubStringMaybeCrunch;

      EnumOutArcs := EnumOutArcs;
      EnumInArcs := EnumInArcs;
      EnumPaths := EnumPaths;
      EnumStrings := EnumStrings;
      EnumStates := EnumStates;

      NPrefs := NPrefs;
      NSuffs := NSuffs;
      NPrefLetters := NPrefLetters;
      NSuffLetters := NSuffLetters;
      EnumPrefs := EnumPrefs;
      EnumSuffs := EnumSuffs;
      FirstPrefix := FirstPrefix;
      FirstSuffix := FirstSuffix;
      FullLabel := FullLabel;
      PrintPrefs := PrintPrefs;
      PrintSuffs := PrintSuffs;
      
      Build := Build;

      NStates := NStates;
      NArcs := NArcs;
      Count := Count;

      MaxAllocState := MaxAllocState;
      Expand := Expand;
      Crunch := Crunch;
    END;

CONST
  NullLetter = FIRST(Basics.Symbol);
  (*
    A transition from "s" through NullLetter to NullState is used
    to indicate that "s" is final. *)

PROCEDURE New(size: POS := 1): T =
  BEGIN
    WITH
      dag = DAG.New(size)
    DO
      RETURN FromDAG(dag, doc := "", root := NullState)
    END
  END New;

PROCEDURE Copy(
    aut: T;
    size: POS := 1;
    VAR (*OUT*) map: REF ARRAY OF State;
  ): T RAISES {Full} =
  BEGIN
    WITH
      newAut = New(size),
      newRoot = CopyStates(
        from := aut,
        to := newAut,
        s := aut.Root(),
        map := (*OUT*) map
      )
    DO
      newAut.SetRoot(newRoot);
      RETURN newAut
    END
  END Copy;

PROCEDURE CopyStates(
    from: T;
    to: T;
    s: State;
    VAR (*IO*) map: REF ARRAY OF State;
  ): State RAISES {Full} =
  BEGIN
    RETURN DAG.Copy(from.dag, to.dag, s, map)
  END CopyStates;

(************************)
(* PRIMITIVE METHODS    *)
(************************)

PROCEDURE HasArcs(<*UNUSED*> aut: T; s: State): BOOL =
  BEGIN
    RETURN s # NullState AND s # UnitState
  END HasArcs;

PROCEDURE Last(aut: T; s: State): Arc =
  BEGIN
    <* ASSERT s # NullState AND s # UnitState *>
    WITH a = aut.dag.Last(s) DO
      RETURN Arc{symbol := a.rd, dest := a.dest}
    END
  END Last;

PROCEDURE Rest(aut: T; s: State): State =
  BEGIN
    <* ASSERT s # NullState AND s # UnitState *>
    RETURN aut.dag.Rest(s)
  END Rest;

PROCEDURE Append(aut:T; s: State; symbol: Symbol; dest: State): State RAISES {Full} =
  BEGIN
    IF s # NullState THEN
      <* ASSERT symbol > aut.dag.Last(s).rd *>
    END;
    IF dest = NullState THEN
      RETURN s
    ELSE
      WITH
        oldMax = aut.MaxState(),
        t = aut.dag.Append(DAG.Arc{rd := symbol, wr := 0, dest := dest}, s)
      DO
        IF t > oldMax THEN
          (* State is truly new, must compute aut.nSuffs[t], aut.nSuffLetters[t] *)
          WITH ns = aut.nSuffs^, nl = aut.nSuffLetters^ DO
            ns[t] := ns[dest] + ns[s];
            nl[t] := nl[dest] + ns[dest] + nl[s];
          END;
        END;
        RETURN t
      END
    END
  END Append;

PROCEDURE Root(aut: T): State =
  BEGIN
    RETURN aut.root
  END Root;

PROCEDURE SetRoot(aut: T; s: State) =
  BEGIN
    aut.root := s
  END SetRoot;

PROCEDURE MaxState(aut: T): POS =
  BEGIN
    RETURN aut.dag.MaxState()
  END MaxState;

PROCEDURE MaxAllocState(aut: T): POS =
  BEGIN
    RETURN aut.dag.MaxAllocState()
  END MaxAllocState;

PROCEDURE Expand(aut: T; newSize: NAT := 0) RAISES {Full} =
  BEGIN
    aut.dag.Expand(newSize);
    ComputeNSuffs(aut);
  END Expand;

PROCEDURE Crunch(aut: T; keep: REF ARRAY OF State := NIL) =
  VAR nKeep: NAT;
  BEGIN
    DiscardPrefixData(aut);
    IF keep = NIL THEN nKeep := 0 ELSE nKeep := NUMBER(keep^) END;
    WITH
      roots = NEW(REF ARRAY OF State, nKeep+1)^
    DO
      roots[0] := aut.root;
      IF keep # NIL THEN SUBARRAY(roots, 1, nKeep) := keep^ END;
      aut.dag.Crunch(roots);
      aut.root := roots[0];
      IF keep # NIL THEN keep^ := SUBARRAY(roots, 1, nKeep) END;
    END;
    (* Recreate unit state, which may have been crunched out: *)
    MakeUnitState(aut.dag);
    (* Recompute suffix counts: *)
    ComputeNSuffs(aut);
  END Crunch;

(******************************)
(* DERIVED METHODS            *)
(******************************)

PROCEDURE Final(aut: T; s: State): BOOL =
  VAR t: State := s;
  BEGIN
    WHILE t # NullState AND t # UnitState DO t := aut.Rest(t) END;
    RETURN t = UnitState
  END Final;

PROCEDURE R(aut: T; s: State; symbol: Symbol): State =
  VAR t: State := s;
  BEGIN
    IF s = NullState THEN RETURN NullState END;
    WHILE t # NullState AND t # UnitState DO
      WITH a = aut.Last(t) DO
        IF a.symbol = symbol THEN
          RETURN a.dest
        ELSIF a.symbol < symbol THEN
          RETURN NullState
        ELSE
          t := aut.Rest(t)
        END
      END
    END;
    RETURN NullState
  END R;

PROCEDURE OutDeg(aut: T; s: State): NAT =
  VAR t: State := s; n: NAT := 0;
  BEGIN
    WHILE t # NullState AND t # UnitState DO INC(n); t := aut.Rest(t) END;
    RETURN n
  END OutDeg;

PROCEDURE InDeg(aut: T; s: State): NAT =
  BEGIN
    IF aut.root = NullState OR s >= aut.root THEN
      RETURN 0
    ELSE
      ComputePrefixData(aut);
      <* ASSERT aut.nPrefs # NIL *>
      RETURN aut.rdag.OutDeg(aut.rev[s])
    END
  END InDeg;

PROCEDURE First(aut: T; s: State): Arc =
  VAR p, t: State;
  BEGIN
    <* ASSERT s # NullState AND s # UnitState *>
    t := s;
    REPEAT p := t; t := aut.Rest(t) UNTIL t = NullState OR t = UnitState;
    RETURN aut.Last(p)
  END First;

PROCEDURE SetArc(aut: T; s: State; symbol: Symbol; dest: State): State RAISES {Full} =

  PROCEDURE DoSet(t: State): State RAISES {Full} =
    BEGIN
      IF t = NullState OR t = UnitState THEN
        RETURN aut.Append(t, symbol, dest)
      ELSE
        WITH a = aut.Last(t) DO
          IF a.symbol < symbol THEN
            RETURN aut.Append(t, symbol, dest)
          ELSIF a.symbol = symbol THEN
            IF a.dest = dest THEN
              RETURN t
            ELSE
              RETURN aut.Append(aut.Rest(t), symbol, dest)
            END;
          ELSE (* a.symbol > symbol *)
            WITH
              rest = aut.Rest(t),
              newRest = DoSet(rest)
            DO
              IF rest = newRest THEN
                RETURN t
              ELSE
                RETURN aut.Append(newRest, a.symbol, a.dest)
              END
            END
          END
        END
      END
    END DoSet;

  BEGIN
    RETURN DoSet(s)
  END SetArc;

PROCEDURE SetFinal(aut: T; s: State; final: BOOL): State RAISES {Full} =

  PROCEDURE DoSet(t: State): State RAISES {Full} =
    BEGIN
      IF t = NullState OR t = UnitState THEN
        IF final THEN 
          RETURN UnitState 
        ELSE
          RETURN NullState
        END
      ELSE
        WITH 
          a = aut.Last(t),
          rest = aut.Rest(t),
          newRest = DoSet(rest)
        DO
          IF newRest = rest THEN 
            RETURN t
          ELSE
            RETURN aut.Append(newRest, a.symbol, a.dest)
          END
        END
      END
    END DoSet;

  BEGIN
    RETURN DoSet(s)
  END SetFinal;

PROCEDURE Walk(aut: T; s: State; READONLY w: String): State =
  VAR t: State := s;
  BEGIN
    FOR i := 0 TO LAST(w) DO
      <* ASSERT w[i] # NullLetter *>
      IF t = NullState THEN
        RETURN NullState
      ELSE
        t := aut.R(t, w[i])
      END
    END;
    RETURN t
  END Walk;

PROCEDURE Accepts(aut: T; s: State; READONLY w: String): BOOL =
  BEGIN
    RETURN aut.Final(aut.Walk(s, w))
  END Accepts;

PROCEDURE Rank(aut: T; s: State; READONLY w: String; reverse: BOOL := FALSE): NAT =
  BEGIN
    IF reverse THEN
      RETURN RankFromBottom(aut, s, w)
    ELSE
      RETURN RankFromTop(aut, s, w)
    END
  END Rank;

PROCEDURE RankFromTop(aut: T; s: State; READONLY w: String): NAT =
(*
  Computes Rank(aut, s, w, reverse := FALSE).
  *)
  VAR rank: NAT := 0;
      i: NAT := 0;
      t: State := s;
  BEGIN
    WHILE t # NullState AND i <= LAST(w) DO
      (* Tally the empty string, if it is a suffix of "t": *)
      IF aut.Final(t) THEN INC(rank) END;
      WITH x = w[i] DO
        (* Tally all suffixes of "t" starting with symbols less than "x": *)
        <* ASSERT x # NullLetter *>
        IF x > FIRST(Symbol) THEN
          rank := rank + AddNSuffsOfChildren(aut, t, lo := FIRST(Symbol), hi := x-1)
        END;
        (* Now add the rank of "w[i+1..]" in Suff(R(t, x)): *)
        t := aut.R(t, x);
        INC(i);
      END
    END;
    RETURN rank
  END RankFromTop;

PROCEDURE RankFromBottom(aut: T; s: State; READONLY w: String): NAT =
(*
  Computes Rank(aut, s, w, reverse := TRUE).
  *)
  VAR rank: NAT := 0;
      i: NAT := 0;
      t: State := s;
  BEGIN
    WHILE t # NullState AND i <= LAST(w) DO
      WITH x = w[i] DO
        (* Tally all suffixes of "t" starting with symbols greater than "x": *)
        <* ASSERT x # NullLetter *>
        IF x < LAST(Symbol) THEN
          rank := rank + AddNSuffsOfChildren(aut, t, lo := x+1, hi := LAST(Symbol))
        END;
        (* Now add the reverse rank of "w[i+1..]" in Suff(R(t, x)): *)
        t := aut.R(t, x);
        INC(i);
      END
    END;
    IF t # NullState THEN
      (* Tally all suffixes of "t" that are strictly grater than (): *)
      rank := rank + aut.NSuffs(t);
      IF aut.Final(t) THEN DEC(rank) END
    END;
    RETURN rank
  END RankFromBottom;

PROCEDURE AddNSuffsOfChildren(aut: T; s: State; lo, hi: Symbol): NAT =
(*
  Returns the sum of NSuffs(R(s, x)) for all "x" in the range [lo..hi]. *)
  VAR sum: NAT := 0;
      t: State := s;
  BEGIN
    WHILE aut.HasArcs(t) DO
      WITH a = aut.Last(t) DO
        IF a.symbol < lo THEN
          RETURN sum
        ELSIF a.symbol <= hi THEN
          INC(sum, aut.NSuffs(a.dest))
        END;
      END;
      t := aut.Rest(t)
    END;
    RETURN sum
  END AddNSuffsOfChildren;

PROCEDURE AddString(aut: T; s: State; READONLY w: String): State RAISES {Full} =
  BEGIN
    RETURN AddSubString(aut, s, w, TRUE)
  END AddString;

PROCEDURE SubString(aut: T; s: State; READONLY w: String): State RAISES {Full} =
  BEGIN
    RETURN AddSubString(aut, s, w, FALSE)
  END SubString;

PROCEDURE AddSubString(
    aut: T;
    s: State;
    READONLY w: String;
    add: BOOL;
  ): State RAISES {Full} =
(*
  Does AddString(aut, s, w) or SubString(aut, s, w),
  depending on "add". *)

  PROCEDURE DoAddSub(i: NAT; t: State): State RAISES {Full} =
  (*
    Returns AddSubString(aut, t, SUBARRAY(w, i), add). *)
    BEGIN
      IF i = NUMBER(w) THEN
        RETURN aut.SetFinal(t, add)
      ELSE
        <* ASSERT w[i] # NullLetter *>
        RETURN aut.SetArc(t, w[i], DoAddSub(i+1, aut.R(t, w[i])))
      END
    END DoAddSub;

  BEGIN
    RETURN DoAddSub(0, s)
  END AddSubString;

PROCEDURE EnumOutArcs(aut: T; s: State; action: ArcAction) RAISES {Abort} =
  VAR i: NAT := 0;

  PROCEDURE DoEnum(r: State) RAISES {Skip, Abort} =
  (*
    Does EnumOutArcs on a given prefix "r" of the arcs out of "s".
    Raises Skip, Abort iff "action" does so. *)
    BEGIN
      IF r = UnitState OR r = NullState THEN
        (* Ok *)
      ELSE
        DoEnum(aut.Rest(r));
        WITH a = aut.Last(r) DO
          action(i, a);
          INC(i)
        END
      END
    END DoEnum;

  BEGIN
    TRY DoEnum(s) EXCEPT Skip => (* Ok *) END;
  END EnumOutArcs;

PROCEDURE AddStringMaybeCrunch(
    aut: T; 
    VAR s: State;
    READONLY w: String; 
    keep: REF States := NIL;
  ): State =
  BEGIN
    RETURN AddSubStringMaybeCrunch(aut, s, w, TRUE, keep)
  END AddStringMaybeCrunch;

PROCEDURE SubStringMaybeCrunch(
    aut: T; 
    VAR s: State; 
    READONLY w: String; 
    keep: REF States := NIL;
  ): State =
  BEGIN
    RETURN AddSubStringMaybeCrunch(aut, s, w, FALSE, keep)
  END SubStringMaybeCrunch;

PROCEDURE AddSubStringMaybeCrunch(
    aut: T; 
    VAR s: State; 
    READONLY w: String; 
    add: BOOL;
    keep: REF States := NIL;
  ): State =

  PROCEDURE CrunchAndExpandIt() =
  (*
    Does aut.Crunch(s & keep), and perhaps aut.Expand(). *)
    VAR nKeep: CARDINAL := 0;
    BEGIN
      IF keep # NIL THEN nKeep := NUMBER(keep^) END;
      WITH 
        tKeep = NEW(REF States, nKeep + 1) 
      DO
        tKeep[0] := s;
        IF keep # NIL THEN SUBARRAY(tKeep^, 1, nKeep) := keep^ END;
        aut.Crunch(tKeep);
        IF keep # NIL THEN keep^ := SUBARRAY(tKeep^, 1, nKeep) END;
        s := tKeep[0];
      END;
      IF aut.MaxState() * 4 >= aut.MaxAllocState() * 3 THEN 
        ExpandIt()
      END;
    END CrunchAndExpandIt;

  PROCEDURE ExpandIt() =
  (*
    Expands the automaton: *)
    BEGIN
      WITH
        oldSize = aut.MaxAllocState() + 1,
        newSize = ComputeNewSize(oldSize)
      DO
        <* FATAL Basics.Full *>
        BEGIN
          aut.Expand(newSize);
        END;
      END;
    END ExpandIt;

  PROCEDURE ComputeNewSize(oldSize: NAT): NAT =
  (*
    Chooses a suitable new size for automatic expansion.
    *)
    BEGIN
      RETURN
        MIN(
          MIN(
                1 + oldSize + oldSize,
            10001 + oldSize + oldSize DIV 2
          ),
          MIN(
            20001 + oldSize + oldSize DIV 3,
            50001 + oldSize + oldSize DIV 4
          )
        )
    END ComputeNewSize;

  BEGIN
    LOOP
      TRY
        RETURN AddSubString(aut, s, w, add);
      EXCEPT 
        Full => CrunchAndExpandIt();
      END;
    END;
  END AddSubStringMaybeCrunch;

PROCEDURE EnumInArcs(aut: T; s: State; action: ArcAction) RAISES {Abort} =
  VAR i: NAT := 0;
      rr: DAG.State;
  BEGIN
    IF s = NullState AND s >= aut.root THEN
      RETURN
    ELSE
      ComputePrefixData(aut);
      rr := aut.rev[s];
      i := 0;
      TRY
        WITH rdag = aut.rdag, dir = aut.dir^ DO
          WHILE rr # DAG.NullState DO
            WITH a = rdag.Last(rr) DO
              action(i, Arc{symbol := a.rd, dest := dir[a.dest]});
              INC(i);
              rr := rdag.Rest(rr)
            END;
          END;
        END
      EXCEPT
      | Skip => (* Ok *)
      END;
    END;
  END EnumInArcs;

PROCEDURE EnumPaths(
    aut: T;
    s: State;
    enter: StateAction := NIL;
    push: PathAction := NIL;
    pop: PathAction := NIL;
    exit: StateAction := NIL;
  ) RAISES {Abort} =

  VAR len: NAT := 0;

  PROCEDURE DoEnumPaths(t: State) RAISES {Abort} =
  (*
    Does EnumPaths starting at "t", a generic sucessor
    of "s". Assumes "t" is not NullState. *)

    VAR final: BOOL; (* Final(t); set by "EnumRest" at the end of the recursion. *)
        i: NAT := 0; (* Arc index (not including the NullLetter, if any) *)

    PROCEDURE EnumRest(r: State) RAISES {Skip, Abort} =
    (*
      Calls "enter" on "t" and enumerates a given prefix "r" of the arcs out of "t".
      Also sets the local variable "final" of DoEnumStates(t)) to Final(t).
      Raises Skip iff "enter" or "pop" raised "Skip". *)

      BEGIN
        IF r = UnitState OR r = NullState THEN
          final := (r = UnitState);
          IF enter # NIL THEN enter(len, t, final) END;
        ELSE
          EnumRest(aut.Rest(r));
          WITH a = aut.Last(r) DO
            TRY
              IF push # NIL THEN push(len, t, i, a) END;
              INC(len);
              DoEnumPaths(a.dest);
              DEC(len);
            EXCEPT
              Skip => (* Ok *)
            END;
            IF pop # NIL THEN pop(len, t, i, a) END;
            INC(i);
          END
        END
      END EnumRest;

    BEGIN
      <* ASSERT t # NullState *>
      TRY EnumRest(t) EXCEPT Skip => (* Ok *) END;
      IF exit # NIL THEN
        <* FATAL Skip *>
        BEGIN
          exit(len, t, final) (* Shouldn't raise "Skip" *)
        END
      END
    END DoEnumPaths;

  BEGIN
    IF s # NullState THEN DoEnumPaths(s) END;
  END EnumPaths;

PROCEDURE EnumStrings(
    aut: T;
    s: State;
    enter: StringAction := NIL;
    exit: StringAction := NIL;
  ) RAISES {Abort} =

  VAR rw: REF String := NEW(REF String, 100);

  PROCEDURE PathEnter (* : StateAction *) (
      len: NAT;
      s: State;
      final: BOOL;
    ) RAISES {Skip, Abort} =
    BEGIN
      IF enter # NIL THEN enter(SUBARRAY(rw^, 0, len), s, final) END
    END PathEnter;

  PROCEDURE PathPush (* : PathAction *) (
      len: NAT;
      <*UNUSED*> org: State;
      <*UNUSED*> i: NAT;
      arc: Arc;
    ) RAISES {} =
    BEGIN
      Basics.ExpandString(rw, len + 1);
      rw[len] := arc.symbol;
    END PathPush;

  PROCEDURE PathPop (* : PathAction *) (
      len: NAT;
      <*UNUSED*> org: State;
      <*UNUSED*> i: NAT;
      arc: Arc;
    ) RAISES {} =
    BEGIN
      <* ASSERT rw[len] = arc.symbol *>
    END PathPop;

  PROCEDURE PathExit (* : StateAction *) (
      len: NAT;
      s: State;
      final: BOOL;
    ) RAISES {Skip, Abort} =
    BEGIN
      IF exit # NIL THEN exit(SUBARRAY(rw^, 0, len), s, final) END
    END PathExit;

  BEGIN
    aut.EnumPaths(
      s,
      enter := PathEnter,
      push := PathPush,
      pop := PathPop,
      exit := PathExit
    );
  END EnumStrings;

PROCEDURE EnumStates(
    aut: T;
    READONLY base: ARRAY OF State;
    substates: BOOL := FALSE;
    enter: StateAction := NIL;
    exit: StateAction := NIL;
  ) RAISES {Abort} =

  VAR maxState: State;
  BEGIN
    (* Computes maximum reachable state: *)
    maxState := NullState;
    FOR i := 0 TO LAST(base) DO maxState := MAX(maxState, base[i]) END;

    WITH
      len = NEW(REF ARRAY OF NAT, maxState + 1)^
    DO

      (* Initialize path lengths: *)
      FOR s := 0 TO maxState DO len[s] := LAST(NAT) END;
      FOR i := 0 TO LAST(base) DO len[base[i]] := 0 END;

      (* First pass: scan states from highest to lowest,
      propagating "len" and "enter"ing all reachable states: *)

      FOR s := maxState TO 1 BY -1 DO
        WITH
          ls = len[s],
          final = aut.Final(s)
        DO
          IF ls < LAST(NAT) THEN
            TRY
              IF enter # NIL THEN enter(ls, s, final) END;
              IF substates THEN
                IF s # UnitState THEN
                  WITH d = aut.Last(s).dest, r = aut.Rest(s) DO
                    len[d] := MIN(len[d], ls + 1);
                    len[r] := MIN(len[r], ls + 1)
                  END
                END
              ELSE
                VAR t := s;
                BEGIN
                  WHILE t # UnitState AND t # NullState DO
                    WITH d = aut.Last(t).dest DO
                      len[d] := MIN(len[d], ls + 1)
                    END;
                    t := aut.Rest(t)
                  END
                END
              END;
            EXCEPT
            | Skip => (* Ignore any arcs out of "s" *)
            END
          END
        END
      END;

      (* Second pass: scan states from low to high, "exit"ing all
      reachable ones: *)

      FOR s := 1 TO maxState DO
        WITH
          ls = len[s],
          final = aut.Final(s)
        DO
          IF ls < LAST(NAT) THEN
            TRY
              IF exit # NIL THEN exit(ls, s, final) END;
            EXCEPT
            | Skip => <* ASSERT FALSE *>
            END;
          END
        END
      END;
    END
  END EnumStates;

CONST
  DumpHeader = "Reduced.Dump (format of 91-12-21)";

PROCEDURE Dump(wr: Wr.T; aut: T) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(wr, "Begin " & DumpHeader); Wr.PutChar(wr, '\n');
    DumpDoc(wr, aut.doc);
    Wr.PutText(wr, "root = " & Fmt.Int(aut.root) & "\n");
    DAG.Dump(wr, aut.dag);
    Wr.PutText(wr, "End " & DumpHeader); Wr.PutChar(wr, '\n');
    Wr.Flush(wr);
  END Dump;
  
EXCEPTION  (* Syntax errors in dump file: *)
  MissingFinalNewLine;
  InvalidHeader;
  InvalidFooter;

PROCEDURE Load(rd: Rd.T; minSize: POS := 1): T =
  <* FATAL Rd.Failure, Rd.EndOfFile, Thread.Alerted *>
  <* FATAL InvalidHeader, InvalidFooter *>
  BEGIN
    WITH hdr = Rd.GetLine(rd) DO
      IF NOT Text.Equal(hdr, "Begin " & DumpHeader) THEN RAISE InvalidHeader END;
    END;
    WITH
      doc = LoadDoc(rd),
      root = ReadParam(rd, "root = "),
      dag = DAG.Load(rd, minSize)
    DO
      WITH hdr = Rd.GetLine(rd) DO
        IF NOT Text.Equal(hdr, "End " & DumpHeader) THEN RAISE InvalidFooter END
      END;
      RETURN FromDAG(dag, doc, root)
    END;
  END Load;
  
PROCEDURE Print(
    wr: Wr.T;
    aut: T;
    e: Encoding.T;
  ) =
  CONST FinalCode = ARRAY BOOL OF CHAR {' ', '*'};

  PROCEDURE PrintState (* : StateAction *) (
      <*UNUSED*> len: NAT;
      s: State;
      final: BOOL;
    ) RAISES {} =
    <* FATAL Wr.Failure, Thread.Alerted *>
    BEGIN
      WrNat(wr, s);
      Wr.PutChar(wr, ' ');
      Wr.PutChar(wr, FinalCode[final]);
      Wr.PutChar(wr, '\n');
      <* FATAL Skip, Abort *>
      BEGIN
        aut.EnumOutArcs(s, PrintArc);
      END
    END PrintState;

  PROCEDURE PrintArc (* : ArcAction *) (
      <*UNUSED*> i: NAT;
      arc: Arc;
    ) RAISES {} =
    <* FATAL Wr.Failure, Thread.Alerted *>
    BEGIN
      Wr.PutChar(wr, ' ');
      Wr.PutChar(wr, ' ');
      e.PrintLetter(wr, arc.symbol);
      Wr.PutChar(wr, ' ');
      Wr.PutChar(wr, '-');
      Wr.PutChar(wr, '>');
      Wr.PutChar(wr, ' ');
      WrNat(wr, arc.dest);
      Wr.PutChar(wr, '\n');
    END PrintArc;

  <* FATAL Wr.Failure, Thread.Alerted, Abort *>
  BEGIN
    Wr.PutText(wr, "root = " & Fmt.Int(aut.root) & "\n");
    Wr.PutChar(wr, '\n');
    aut.EnumStates(base := ARRAY OF State{aut.root}, enter := PrintState);
    Wr.Flush(wr);
  END Print;

PROCEDURE NStates(aut: T; READONLY base: ARRAY OF State): NAT =
  VAR nStates: NAT := 0;

  PROCEDURE CountState (* : StateAction *) (
      <*UNUSED*> len: NAT;
      <*UNUSED*> s: State;
      <*UNUSED*> final: BOOL;
    ) RAISES {} =
    BEGIN
      INC(nStates);
    END CountState;

  <* FATAL Abort *>
  BEGIN
    aut.EnumStates(base := base, enter := CountState);
    RETURN nStates;
  END NStates;

PROCEDURE NArcs(aut: T; READONLY base: ARRAY OF State): NAT =
  VAR nArcs: NAT := 0;

  PROCEDURE CountArcs (* : StateAction *) (
      <*UNUSED*> len: NAT;
      s: State;
      <*UNUSED*> final: BOOL;
    ) RAISES {} =
    BEGIN
      INC(nArcs, aut.OutDeg(s));
    END CountArcs;

  <* FATAL Abort *>
  BEGIN
    aut.EnumStates(base := base, enter := CountArcs);
    RETURN nArcs;
  END NArcs;

PROCEDURE Count(aut: T; READONLY base: ARRAY OF State): Counts =

  VAR nStates: NAT := 0;
      nArcs: NAT := 0;
      nFinals: NAT := 0;
      nSubStates: NAT := 0;
      nStrings: NAT := 0;
      nLetters: NAT := 0;

  PROCEDURE CountState (* : StateAction *) (
      <*UNUSED*> len: NAT;
      s: State;
      final: BOOL;
    ) RAISES {} =
    BEGIN
      INC(nStates);
      INC(nArcs, aut.OutDeg(s));
      IF final THEN INC(nFinals) END;
    END CountState;

  PROCEDURE CountSubState (* : StateAction *) (
      <*UNUSED*> len: NAT;
      <*UNUSED*> s: State;
      <*UNUSED*> final: BOOL;
    ) RAISES {} =
    BEGIN
      INC(nSubStates);
    END CountSubState;

  <* FATAL Abort *>
  BEGIN
    (*
      Note that "nStates" counts only the DAG states reachable by "R"
      chains, not those reachable by "Last-Rest" chains.

      Note also that "nArcs" is the sum of the outdegre of all
      "R"-reachable states of the Reduced.T; thus, DAG arcs that
      are shared by more than one reachable state are counted
      more than once.
    *)
    aut.EnumStates(base := base, enter := CountState);
    aut.EnumStates(base := base, substates := TRUE, enter := CountSubState);
    nStrings := 0;
    nLetters := 0;
    FOR i := 0 TO LAST(base) DO
      nStrings := nStrings + aut.NSuffs(base[i]);
      nLetters := nLetters + aut.NSuffLetters(base[i]);
    END;
    RETURN Counts{
      strings := nStrings,
      symbols := nLetters,
      states := nStates,
      substates := nSubStates,
      arcs := nArcs,
      finals := nFinals
    }
  END Count;

PROCEDURE NPrefs(aut: T; s: State): NAT =
  BEGIN
    IF s = NullState OR s > aut.root THEN
      RETURN 0
    ELSIF s = aut.root THEN
      RETURN 1
    ELSE
      ComputePrefixData(aut);
      <* ASSERT aut.nPrefs # NIL *>
      WITH np = aut.nPrefs^ DO
        <* ASSERT LAST(np) >= aut.root *>
        RETURN np[s]
      END;
    END
  END NPrefs;

PROCEDURE NSuffs(aut: T; s: State): NAT =
  BEGIN
    IF s = NullState THEN
      RETURN 0
    ELSIF s = UnitState THEN
      RETURN 1
    ELSE
      <* ASSERT aut.nSuffs # NIL *>
      WITH ns = aut.nSuffs^ DO
        <* ASSERT LAST(ns) >= s *>
        RETURN ns[s]
      END;
    END;
  END NSuffs;

PROCEDURE NPrefLetters(aut: T; s: State): NAT =
  BEGIN
    IF s = NullState OR s > aut.root THEN
      (* No paths from "aut.root" to "s":  *)
      RETURN 0
    ELSIF s = aut.root THEN
      (* One path from "aut.root" to "s", of length 0: *)
      RETURN 0
    ELSE
      ComputePrefixData(aut);
      <* ASSERT aut.nPrefLetters # NIL *>
      WITH nl = aut.nPrefLetters^ DO
        <* ASSERT LAST(nl) >= aut.root *>
        RETURN nl[s]
      END;
    END
  END NPrefLetters;

PROCEDURE NSuffLetters(aut: T; s: State): NAT =
  BEGIN
    IF s = NullState OR s = UnitState THEN
      RETURN 0
    ELSE
      <* ASSERT aut.nSuffLetters # NIL *>
      WITH nl = aut.nSuffLetters^ DO
        <* ASSERT LAST(nl) >= s *>
        RETURN nl[s]
      END;
    END;
  END NSuffLetters;

PROCEDURE EnumPrefs(aut: T; s: State; action: PrefixAction) RAISES {Abort} =
  VAR rs: REF String := NEW(REF String, 100);

  PROCEDURE PrefixPush (* : DAG.PathAction *) (
      len: NAT;
      <*UNUSED*> org: State;
      <*UNUSED*> i: NAT;
      a: DAG.Arc
    ) RAISES {Skip, Abort} =
    BEGIN
      IF a.rd = NullLetter THEN
        action(SUBARRAY(rs^, 0, len));
        RAISE Skip
      ELSE
        Basics.ExpandString(rs, len+1);
        rs[len] := a.rd
      END
    END PrefixPush;

  BEGIN
    IF s = NullState OR s > aut.root THEN
      RETURN
    ELSIF s = aut.root THEN
      action(String{});
    ELSE
      ComputePrefixData(aut);
      aut.rdag.EnumPaths(aut.rev[s], push := PrefixPush)
    END;
  END EnumPrefs;

PROCEDURE EnumSuffs(aut: T; s: State; action: SuffixAction) RAISES {Abort} =

  PROCEDURE SuffixAction (* : StringAction *) (
       READONLY w: String;
       <*UNUSED*> dest: State;
       final: BOOL;
    ) RAISES {Abort} =
    BEGIN
      IF final THEN action(w) END
    END SuffixAction;

  BEGIN
    aut.EnumStrings(s, enter := SuffixAction)
  END EnumSuffs;

PROCEDURE PrintPrefs(
    aut: T;
    s: State;
    spr: StringPrinter.T;
  ) =

  PROCEDURE PrintPrefix (* : PrefixAction *) (
      READONLY w: String
    ) RAISES {Abort} =
    BEGIN
      spr.PutString(w, rev := TRUE)
    END PrintPrefix;

  BEGIN
    TRY aut.EnumPrefs(s, action := PrintPrefix) EXCEPT Abort => (* Ok *) END;
    spr.Reset();
  END PrintPrefs;

PROCEDURE PrintSuffs(
    aut: T;
    s: State;
    spr: StringPrinter.T;
  ) =

  PROCEDURE PrintSuffix (* : SuffixAction *) (
      READONLY w: String
    ) RAISES {Abort} =
    BEGIN
      spr.PutString(w, rev := FALSE)
    END PrintSuffix;

  BEGIN
    TRY aut.EnumSuffs(s, action := PrintSuffix) EXCEPT Abort => (* Ok *) END;
    spr.Reset();
  END PrintSuffs;

PROCEDURE FirstPrefix(aut: T; s: State; e: Encoding.T): TEXT =

  VAR wr := TextWr.New();

  PROCEDURE FP(t: DAG.State) =
  (*
    Prints into "wr" the string of the first path from "t"
    to an initial state (recognized by its first arc being labelled
    with NullLetter). *)
    <* FATAL Wr.Failure, Thread.Alerted *>
    BEGIN
      <* ASSERT t # NullState *>
      WITH a = aut.rdag.First(t) DO
        IF a.rd = NullLetter THEN
          RETURN
        ELSE
          FP(a.dest);
          e.PrintLetter(wr, a.rd);
        END;
      END;
    END FP;

  BEGIN
    <* ASSERT s # NullState *>
    <* ASSERT s <= aut.root *>
    ComputePrefixData(aut);
    FP(aut.rev[s]);
    RETURN TextWr.ToText(wr)
  END FirstPrefix;

PROCEDURE FirstSuffix(aut: T; s: State; e: Encoding.T): TEXT =
  VAR t: State := s;
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    WITH wr = TextWr.New() DO
      LOOP
        <* ASSERT t # NullState *>
        WITH a = aut.dag.First(t) DO
          IF a.rd = NullLetter THEN RETURN TextWr.ToText(wr) END;
          e.PrintLetter(wr, a.rd);
          t := a.dest
        END
      END;
    END;
  END FirstSuffix;

PROCEDURE FullLabel(aut: T; s: State; e: Encoding.T; sep: TEXT := ":"): TEXT =
  BEGIN
    RETURN
      aut.FirstPrefix(s, e) & sep & aut.FirstSuffix(s, e)
  END FullLabel;
  
VAR (*CONST*) DefaultEncoding: Encoding.T := PlainEncoding.New();

PROCEDURE Build(
    aut: T;
    next: NextStringProc;        (* Client input procedure *)
    wr: Wr.T := NIL;             (* Writer for progress report, etc: *)
    e: Encoding.T := NIL;        (* Symbol/CHAR encoding for printout *)
    reportInterval: POS := 1000; (* Print a report every this many input strings *)
    flagRedundant: BOOL := TRUE; (* TRUE to print warnings on redundant operations *)
  ) RAISES {Abort} =

  VAR
    nStrings: ARRAY BOOL OF NAT := ARRAY OF NAT {0, 0};
    nLetters: ARRAY BOOL OF NAT := ARRAY OF NAT {0, 0};

  VAR
    lastCrunchMaxState: State := NullState;  (* aut.MaxState() after last GC run *)

  PROCEDURE CrunchIt() =
  (*
    Does aut.Crunch(), printing status reports. *)
    <* FATAL Wr.Failure, Thread.Alerted *>
    BEGIN
      WITH size = aut.MaxAllocState() + 1 DO
        IF wr # NIL THEN
          IF NOT TimeToReport() THEN PrintStatusReport() END;
          Wr.PutText(wr, "    * (crunching, alloc = " &  Fmt.Int(size) & "...");
          Wr.Flush(wr);
        END;
        aut.Crunch(pSt);
        IF wr # NIL THEN 
          Wr.PutText(wr, ")\n");
          PrintStatusReport();
        END;
      END;
      lastCrunchMaxState := aut.MaxState();
    END CrunchIt;

  PROCEDURE ExpandIt() =
  (*
    Expands the automaton, printing some noise: *)
    <* FATAL Wr.Failure, Thread.Alerted *>
    BEGIN
      WITH
        oldSize = aut.MaxAllocState() + 1,
        newSize = ComputeNewSize(oldSize)
      DO
        IF wr # NIL THEN
          IF NOT TimeToReport() THEN PrintStatusReport() END;
          Wr.PutText(wr, "    * (expanding from ");
          Wr.PutText(wr, Fmt.Int(oldSize));
          Wr.PutText(wr, " to ");
          Wr.PutText(wr, Fmt.Int(newSize));
          Wr.PutText(wr, "...");
          Wr.Flush(wr);
        END;
        <* FATAL Basics.Full *>
        BEGIN
          aut.Expand(newSize);
        END;
        IF wr # NIL THEN 
          Wr.PutText(wr, ")\n");
          PrintStatusReport();
        END;
      END;
    END ExpandIt;

  PROCEDURE ComputeNewSize(oldSize: NAT): NAT =
  (*
    Chooses a suitable new size for automatic expansion.
    *)
    BEGIN
      RETURN
        MIN(
          MIN(
                1 + oldSize + oldSize,
            10001 + oldSize + oldSize DIV 2
          ),
          MIN(
            20001 + oldSize + oldSize DIV 3,
            50001 + oldSize + oldSize DIV 4
          )
        )
    END ComputeNewSize;

  VAR
    pN: NAT := 0;                              (* Number of pending SetArc actions *)
    pSt: REF States := NEW(REF States, 100);   (* State arguments for pending actions *)
    pLet: REF String := NEW(REF String, 100);  (* Symbol arguments for pending SetArcs *)
      (*
        For efficiency reasons, the strings returned by "next" are not added or
        deleted right away.  Instead, we keep a list "pa[0..pN-1]" of "pending
        SetArc actions" to be performed later.
        
        The "i"th pending SetArc action is described by the pair /pSt[i],
        pLet[i]". The action consists in computing a new state "new[i]/ that 
        is like "pSt[i]" except that under the symbol
        "pLet[i]" it goes to the state "new[i+1]". As a special case,
        "new[pN]" is defined to be "pSt[pN]" itself, and "pLet[pN]" is
        not used.
        
        Thus, to flush out all pending actions we must repeat
        
|          pSt[pN-1] := aut.SetArc(pSt[pN-1], pLet[pN-1], s);
|          DEC(pN)

        until pN = 0, and finally set the automaton's root to the
        resuting state "pSt[0]". *)
    
  PROCEDURE DoSetFinal(s: State; final: BOOL): State =
  (*
    Computes aut.SetFinal(s, final), but expands "aut" if necessary
    (instead of raising "Full"). *)
    BEGIN
      LOOP
        TRY
          RETURN aut.SetFinal(s, final);
        EXCEPT 
          Full => ExpandIt();
        END;
      END;
    END DoSetFinal;

  PROCEDURE DoOnePendingAction() =
  (*
    Performs the last pending "SetArc" action. Expands if necessary. 
    *)
    BEGIN
      <* ASSERT pN > 0 *>
      WITH 
        st = pSt^, let = pLet^, 
        n1 = pN-1 
      DO
        LOOP
          TRY
            st[n1] := aut.SetArc(st[n1], let[n1], st[pN]);
            pN := n1;
            RETURN
          EXCEPT 
            Full => ExpandIt();
          END;
        END;
      END
    END DoOnePendingAction;

  PROCEDURE FlagRedundant(READONLY w: String; add: BOOL) =
    <* FATAL Wr.Failure, Thread.Alerted *>
    BEGIN
      Wr.PutText(wr, "  ** redundant command: ");
      PrintCommand(w, add);
      Wr.PutChar(wr, '\n');
    END FlagRedundant;

  PROCEDURE AddOrSubIt(READONLY w: String; add: BOOL := TRUE) =
  (*
    Adds or deletes "w", crunching and/or expanding if necessary:
    *)
    BEGIN
      (* Crunch the automaton, if it looks worth it: *)
      WITH alloc = aut.MaxAllocState() DO
        IF aut.MaxState() + 30*NUMBER(w) >= alloc
        AND (aut.MaxState() - lastCrunchMaxState) >= alloc DIV 10
        THEN
          (* We are close to the allocated size,
            and about 10% of the current dag states were
            created since the last CG run. Better run GC again... *)
          CrunchIt();
        END;
      END;

      (* Now add the string, and expand if doesn't fit: *)
      VAR np: CARDINAL := 0;
      BEGIN
        (* Skip pending actions that match the symbols of w: *)
        WITH maxp = MIN(pN, NUMBER(w)), let = pLet^ DO
          WHILE np < maxp AND w[np] = let[np] DO INC(np) END
        END;
        
        (* Flush any remaining actions: *)
        IF pN > np THEN
          REPEAT DoOnePendingAction() UNTIL pN <= np
        END;
        
        (* Stack new actions corresponding to remaining symbols of w: *)
        Basics.ExpandString(pLet, NUMBER(w));
        ExpandStates(pSt, NUMBER(w) + 1);
        <* ASSERT pN = np *>
        WITH st = pSt^, let = pLet^ DO
          FOR i := np TO LAST(w) DO
            let[i] := w[i];
            st[i+1] := aut.R(st[i], w[i]);
          END;
          pN := NUMBER(w);
          IF add = aut.Final(pSt[pN]) THEN
            (* Command was superflous *)
            IF wr # NIL AND flagRedundant THEN
              FlagRedundant(w, add)
            END;
          ELSE
            pSt[pN] := DoSetFinal(pSt[pN], add)
          END;
        END;
      END;
    END AddOrSubIt;

  VAR
    refw: REF String := NEW(REF String, 100);  (* String buffer *)
    len: NAT;                                  (* Length of string *)
    add: BOOL;                                 (* Add/delete flag for refw^ *)
    
  PROCEDURE PrintStatusReport() =
    <* FATAL Wr.Failure, Thread.Alerted *>
    BEGIN
      Wr.PutText(wr, "    * ");
      WITH 
        ns = nStrings[FALSE] + nStrings[TRUE],
        nl = nLetters[FALSE] + nLetters[TRUE]
      DO
        Wr.PutText(wr, Fmt.Pad(Fmt.Int(ns), 8) &  " strings ");
        Wr.PutText(wr, Fmt.Pad(Fmt.Int(nl), 8) &  " symbols ");
      END;
      Wr.PutText(wr, Fmt.Pad(Fmt.Int(aut.MaxState()), 8) & " dag states  ");
      PrintCommand(SUBARRAY(refw^, 0, len), add);
      Wr.PutText(wr, "\n");
      Wr.Flush(wr);
    END PrintStatusReport;

  PROCEDURE PrintCommand(READONLY w: String; add: BOOL) =
    <* FATAL Wr.Failure, Thread.Alerted *>
    BEGIN
      Wr.PutChar(wr, ARRAY BOOL OF CHAR{'-', '+'}[add]);
      Wr.PutChar(wr, ' ');
      TRY e.PrintString(wr, w) EXCEPT Encoding.BadString => (*IGNORE*) END;
    END PrintCommand;

  PROCEDURE PrintFinalReport() =
    <* FATAL Wr.Failure, Thread.Alerted *>
    BEGIN
      Wr.PutText(wr, "\n");

      Wr.PutText(wr, "Input statistics:\n");
      Wr.PutText(wr, "\n");

      Wr.PutText(wr, "add:  ");
      Wr.PutText(wr, Fmt.Pad(Fmt.Int(nStrings[TRUE]), 8) & " strings ");
      Wr.PutText(wr, Fmt.Pad(Fmt.Int(nLetters[TRUE]), 8) & " symbols ");
      Wr.PutText(wr, "\n");

      Wr.PutText(wr, "sub:  ");
      Wr.PutText(wr, Fmt.Pad(Fmt.Int(nStrings[FALSE]), 8) & " strings ");
      Wr.PutText(wr, Fmt.Pad(Fmt.Int(nLetters[FALSE]), 8) & " symbols ");
      Wr.PutText(wr, "\n");

      Wr.PutText(wr, "------");
      Wr.PutText(wr, "--------");
      Wr.PutText(wr, "---------");
      Wr.PutText(wr, "--------");
      Wr.PutText(wr, "---------");
      Wr.PutText(wr, "\n");
      WITH 
        ns = nStrings[FALSE] + nStrings[TRUE],
        nl = nLetters[FALSE] + nLetters[TRUE]
      DO
        Wr.PutText(wr, "tot:  ");
        Wr.PutText(wr, Fmt.Pad(Fmt.Int(ns), 8) & " strings ");
        Wr.PutText(wr, Fmt.Pad(Fmt.Int(nl), 8) & " symbols ");
      END;
      Wr.PutText(wr, "\n");
      Wr.PutText(wr, "\n");
      WITH ct = aut.Count(ARRAY OF State{aut.Root()}) DO
        PrintCounts(wr, ct)
      END;
    END PrintFinalReport;

  PROCEDURE TimeToReport(): BOOL =
    BEGIN
      RETURN (nStrings[FALSE] + nStrings[TRUE]) MOD reportInterval = 0
    END TimeToReport;
  
  BEGIN (*Build*)
    pN := 0;
    pSt[0] := aut.Root();
    IF wr # NIL THEN
      IF e = NIL THEN e := DefaultEncoding END;
    END;
    TRY
      TRY
        LOOP
          next((*IO*) refw, (*OUT*) len, (*OUT*) add);
          WITH w = SUBARRAY(refw^, 0, len) DO
            AddOrSubIt(w, add);
            INC(nStrings[add]);
            INC(nLetters[add], NUMBER(w));
            IF wr # NIL AND TimeToReport() THEN 
              PrintStatusReport(); 
            END;
          END
        END
      EXCEPT
      | Done => (* Ok *)
      END;
    FINALLY (* Normally, or in case of "Abort" *)
      WHILE pN > 0 DO DoOnePendingAction() END;
      aut.SetRoot(pSt[0])
    END;
    CrunchIt();
    IF wr # NIL THEN PrintFinalReport() END;
  END Build;

(************************)
(* AUXILIARY OPERATIONS *)
(************************)

PROCEDURE FromDAG(dag: DAG.T; doc: TEXT; root: DAG.State): T =
(*
  Bulds a Reduced.T given the underlying DAG and the root state.
  The DAG must satisfy the implementation conventions about the
  use of NullLetter and the ordering of arc labels. *)
  BEGIN
    WITH
      aut = NEW(T,
        doc := doc,
        dag := dag,
        root := root,
        nSuffs := NIL,
        nSuffLetters := NIL,
        prefRoot := NullState,
        rdag := NIL,
        rev := NIL,
        dir := NIL,
        nPrefs := NIL,
        nPrefLetters := NIL
      )
    DO
      MakeUnitState(dag);
      ComputeNSuffs(aut);
      RETURN aut;
    END
  END FromDAG;

PROCEDURE MakeUnitState(dag: DAG.T) =
(*
  Creates the "unit" state in the "dag", if necesary: *)
  <* FATAL Basics.Full *>
  BEGIN
    WITH 
      unit = dag.Append(
        last := DAG.Arc{rd := 0, wr := 0, dest := DAG.NullState},
        rest := DAG.NullState
      )
    DO
      <* ASSERT unit = UnitState *>
    END
  END MakeUnitState;

PROCEDURE DiscardPrefixData(aut: T) =
  BEGIN
    aut.prefRoot := NullState;
    aut.rdag := NIL;
    aut.rev := NIL;
    aut.dir := NIL;
    aut.nPrefs := NIL;
    aut.nPrefLetters := NIL;
  END DiscardPrefixData;

PROCEDURE ComputePrefixData(aut: T) =
(*
  If the prefix tables of "aut" ("rdag", "rev", "dir" and "nPrefs")
  are missing or out of date (e.g., because the root has changed),
  recomputes them for the current root.
  *)
  BEGIN
    IF aut.root = NullState OR aut.prefRoot = aut.root THEN
      (* Current prefix data is still OK *)
      RETURN
    ELSE
      DiscardPrefixData(aut);
      ComputeReverseDAG(aut);
      ComputeNPrefs(aut);
      aut.prefRoot := aut.root;
    END
  END ComputePrefixData;

PROCEDURE ComputeReverseDAG(aut: T) =
(*
  Computes "aut.rdag", "aut.rev", "aut.dir" from "aut.dag" and "aut.root".
  *)
  VAR maxRev: DAG.State;
  <* FATAL Basics.Full *>
  BEGIN
    <* ASSERT aut.root # NullState *>
    (*
      Note that the reverse DAG may be bigger than "aut.dag",
      because one arc of "aut" that is shared by two
      reachable states will give rise to two distinct reverse arcs.
      *)
    WITH
      totArcs = aut.Count(ARRAY OF State{aut.root}).arcs,
      rdag = DAG.New(size := totArcs + 1),
      rev = NEW(REF ARRAY OF DAG.State, aut.root+1),
      r = rev^
    DO
      FOR t := 0 TO aut.root DO r[t] := NullState END;
      (* Create one transition from the root state to NullState with NullLetter: *)
      r[aut.root] := rdag.Append(
        last := DAG.Arc{rd := NullLetter, wr := NullLetter, dest := NullState},
        rest := NullState
      );

      (*
        For every state "s" reachable from "root", and every
        proper successor "t" of "s" in "aut", add to "rdag" an arc from
        "rev[t]" to "rev[s]", and update rev[t] accordingly.
        Note that it is important to process the states from high to low,
        so that rev[s] is stable by the time we add the reverse arcs into it. *)

      maxRev := NullState;
      FOR s := aut.root TO 1 BY -1 DO
        IF r[s] # NullState THEN
          (* "s" is reachable from the root state. *)

          (* Update "maxRev": *)
          maxRev := MAX(maxRev, r[s]);

          (* Enumerate its outgoing arcs, and add them to the reverse dag: *)
          VAR t: State := s;
          BEGIN
            WHILE aut.HasArcs(t) DO
              WITH a = aut.Last(t) DO
                r[a.dest] := rdag.Append(
                  last := DAG.Arc{rd := a.symbol, wr := NullLetter, dest := r[s]},
                  rest := r[a.dest]
                );
                t := aut.Rest(t)
              END;
            END;
          END;
        END;
      END;

      (* Store results in "aut": *)
      aut.rdag := rdag;
      aut.rev := rev;

      (* Build table from reversed state to direct state: *)
      WITH
        dir = NEW(REF ARRAY OF State, maxRev + 1),
        d = dir^
      DO
        FOR t := 0 TO maxRev DO d[t] := NullState END;
        FOR s := 0 TO aut.root DO
          WITH t = r[s] DO
            IF t # NullState THEN
              <* ASSERT d[t] = NullState *>
              d[t] := s
            END
          END
        END;
        aut.dir := dir
      END;
    END;
  END ComputeReverseDAG;

PROCEDURE ComputeNPrefs(aut: T) =
(*
  Computes the tables "nPrefs" and "nPrefLetters", for the current root state,
  from 
  *)
  VAR t: State;
  BEGIN
    <* ASSERT aut.root # NullState *>
    (* Allocates the vector if necessary. *)
    WITH
      minSize = aut.root + 1
    DO
      IF aut.nPrefs = NIL OR NUMBER(aut.nPrefs^) < minSize THEN
        aut.nPrefs := NEW(REF ARRAY OF NAT, minSize)
      END;
      IF aut.nPrefLetters = NIL OR NUMBER(aut.nPrefLetters^) < minSize THEN
        aut.nPrefLetters := NEW(REF ARRAY OF NAT, minSize)
      END;
    END;

    WITH
      np = aut.nPrefs^,
      nl = aut.nPrefLetters^,
      maxState = aut.root
    DO
      FOR s := NullState TO maxState DO 
        np[s] := 0; nl[s] := 0 
      END;
      IF maxState > NullState THEN np[maxState] := 1 END;
      FOR s := maxState TO 1 BY -1 DO
        WITH  nps = np[s], nls = nl[s] DO
          IF nps # 0 THEN
            t := s;
            WHILE aut.HasArcs(t)  DO
              WITH
                d = (aut.Last(t)).dest
              DO
                <* ASSERT d < s *>
                INC(np[d], nps);
                INC(nl[d], nls + nps);
              END;
              t := aut.Rest(t)
            END
          END
        END;
      END
    END
  END ComputeNPrefs;

PROCEDURE ComputeNSuffs(aut: T) =
(*
  Computes the tables "nSuffs" and "nSuffLetters", from "aut.dag". *)
  BEGIN
    (* Allocates the vectors if necessary. *)
    WITH
      minSize = aut.MaxAllocState() + 1
    DO
      IF aut.nSuffs = NIL OR NUMBER(aut.nSuffs^) < minSize THEN
        aut.nSuffs := NEW(REF ARRAY OF NAT, minSize)
      END;
      IF aut.nSuffLetters = NIL OR NUMBER(aut.nSuffLetters^) < minSize THEN
        aut.nSuffLetters := NEW(REF ARRAY OF NAT, minSize)
      END;
    END;

    (* Computes suffix counts and sizes, by a postorder scan *)
    WITH
      ns = aut.nSuffs^,
      nl = aut.nSuffLetters^,
      maxState = aut.MaxState()
    DO
      <* ASSERT maxState >= UnitState *>
      ns[NullState] := 0;  nl[NullState] := 0;
      ns[UnitState] := 1;  nl[UnitState] := 0;
      FOR s := 2 TO maxState DO
        WITH
          d  = (aut.Last(s)).dest,
          r  = aut.Rest(s)
        DO
          <* ASSERT (d < s) AND (r < s) *>
          ns[s] := ns[d] + ns[r];
          nl[s] := nl[d] + ns[d] + nl[r]
        END;
      END;
      (* Clear out remaining counts, just for tidiness: *)
      FOR i := maxState + 1 TO LAST(ns) DO ns[i] := 0; nl[i] := 0 END;
    END;
  END ComputeNSuffs;

CONST DocPrefix: CHAR = '|';
  
PROCEDURE DumpDoc(wr: Wr.T; doc: TEXT) =
(* 
  Writes the given "doc" text to "wr", with a DocPrefix
  in front of every line.  Supplies a final '\n' if the text is 
  non-empty but does not end with newline. *)
  
  VAR rd: Rd.T := TextRd.New(doc);
  
  PROCEDURE CopyLine() RAISES {Rd.EndOfFile} =
  (*
    Copy one line from "rd" to "wr", prefixed by DocPrefix. 
    Supplies a final '\n' if next line exists but does not end with newline.
    Raises Rd.EndOfFile if there are no more lines in "rd". *)
    
    <* FATAL Rd.Failure, Wr.Failure, Thread.Alerted *>
    VAR c: CHAR;
    BEGIN
      c := Rd.GetChar(rd); (* If EOF here, propagate to caller *)
      Wr.PutChar(wr, DocPrefix);
      Wr.PutChar(wr, c);
      WHILE c # '\n' DO
        TRY c := Rd.GetChar(rd) EXCEPT Rd.EndOfFile => c := '\n' END;
        Wr.PutChar(wr, c)
      END
    END CopyLine;

  BEGIN
    TRY LOOP CopyLine() END EXCEPT Rd.EndOfFile => (* Ok *) END;
  END DumpDoc;

PROCEDURE LoadDoc(rd: Rd.T): TEXT =
(*
  Reads zero or more lines from "rd" that begin with DocPrefix, strips the
  leading DocPrefix of each line, and returns the concatenation of those lines
  as a single TEXT with embedded and terminating newline chars. *)

  VAR wr: Wr.T := TextWr.New();

  PROCEDURE CopyLine() RAISES {Rd.EndOfFile} =
  (*
    Copy one DocPrefix line from "rd" to "wr", removing the DocPrefix
    but leaving the final (mandatory) newline.
    Raises Rd.EndOfFile if "rd" is exhausted or the next char is 
    not DocPrefix. *)
    <* FATAL Rd.Failure, Wr.Failure, Thread.Alerted *>
    <* FATAL MissingFinalNewLine *>
    VAR c: CHAR;
    BEGIN
      c := Rd.GetChar(rd); (* If EOF here, propagate to caller *)
      IF c # DocPrefix THEN Rd.UnGetChar(rd); RAISE Rd.EndOfFile END;
      REPEAT
        TRY c := Rd.GetChar(rd) EXCEPT Rd.EndOfFile => RAISE MissingFinalNewLine END;
        Wr.PutChar(wr, c)
      UNTIL c = '\n'
    END CopyLine;

  VAR twr: TextWr.T;
      txt: TEXT;
  BEGIN
    TRY LOOP CopyLine() END EXCEPT Rd.EndOfFile => (* Ok *) END;
    twr := wr;
    txt := TextWr.ToText(twr);
    RETURN txt
  END LoadDoc;  

PROCEDURE PrintCounts(wr: Wr.T; READONLY ct: Counts) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(wr, " strings  symbols    states   finals     arcs  sub-sts lets/arc\n");
    Wr.PutText(wr, "-------- --------  -------- -------- -------- -------- --------\n");
    Wr.PutText(wr, Fmt.Pad(Fmt.Int(ct.strings),   8));
    Wr.PutChar(wr, ' ');
    Wr.PutText(wr, Fmt.Pad(Fmt.Int(ct.symbols),   8));
    Wr.PutChar(wr, ' ');
    Wr.PutChar(wr, ' ');
    Wr.PutText(wr, Fmt.Pad(Fmt.Int(ct.states),    8));
    Wr.PutChar(wr, ' ');
    Wr.PutText(wr, Fmt.Pad(Fmt.Int(ct.finals),    8));
    Wr.PutChar(wr, ' ');
    Wr.PutText(wr, Fmt.Pad(Fmt.Int(ct.arcs),      8));
    Wr.PutChar(wr, ' ');
    Wr.PutText(wr, Fmt.Pad(Fmt.Int(ct.substates), 8));
    Wr.PutChar(wr, ' ');
    Wr.PutText(wr, 
      Fmt.Pad(Fmt.Real(FLOAT(ct.symbols)/FLOAT(ct.arcs), Fmt.Style.Fix, 3), 8)
    );
    Wr.PutChar(wr, '\n');
  END PrintCounts;
  
PROCEDURE ExpandStates(VAR s: REF States; n: CARDINAL) =
  BEGIN
    WITH nold = NUMBER(s^) DO
      IF nold < n THEN
        WITH
          r = NEW(REF States, MAX(n + 10, 2*nold))
        DO
          SUBARRAY(r^, 0, nold) := s^;
          s := r
        END
      END
    END
  END ExpandStates;

BEGIN
  <* ASSERT NullState = DAG.NullState *>
  <* ASSERT NullLetter < FIRST(Symbol) *>
END Reduced.

(***********************************************************************)
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
 
