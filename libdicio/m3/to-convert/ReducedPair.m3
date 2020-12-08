MODULE ReducedPair;

(* See the copyright and disclaimer note at the end of this file. *)

IMPORT Basics, StringPrinter, Reduced;
FROM Basics IMPORT NAT, BOOL, Skip, Abort;
FROM Reduced IMPORT Symbol, String;

REVEAL
  T = Public BRANDED OBJECT
    OVERRIDES

      HasArcs := HasArcs;
      Last := Last;
      Rest := Rest;
      Root := Root;
      Final := Final;

      OutDeg := OutDeg;
      First := First;
      R := R;
      Walk := Walk;

      EnumArcs := EnumArcs;
      EnumPaths := EnumPaths;
      EnumStrings := EnumStrings;

      Accepts := Accepts;
      EnumSuffs := EnumSuffs;
      NSuffs := NSuffs;
      PrintSuffs := PrintSuffs;

    END;

CONST
  RedNull = Reduced.NullState;
  RedUnit = Reduced.UnitState;

PROCEDURE New(a0, a1: Reduced.T): T =
  BEGIN
    RETURN NEW(T, aut := ARRAY Which OF Reduced.T{a0, a1})
  END New;

(* PRIMITIVE METHODS *)

PROCEDURE HasArcs(<*UNUSED*> p: T; s: State): BOOL =
  BEGIN
    RETURN (s[0] # RedNull AND s[0] # RedUnit)
        OR (s[1] # RedNull AND s[1] # RedUnit)
  END HasArcs;

PROCEDURE Last(p: T; s: State): Arc =
  BEGIN
    IF s[0] = RedNull OR s[0] = RedUnit THEN
      IF s[1] = RedNull OR s[1] = RedUnit THEN
        <* ASSERT FALSE *>
      ELSE
        WITH a1 = p.aut[1].Last(s[1]) DO
          RETURN Arc{symbol := a1.symbol, dest := State{RedNull, a1.dest}}
        END
      END;
    ELSE
      IF s[1] = RedNull OR s[1] = RedUnit THEN
        WITH a0 = p.aut[0].Last(s[0]) DO
          RETURN Arc{symbol := a0.symbol, dest := State{a0.dest, RedNull}}
        END
      ELSE
        WITH
          a0 = p.aut[0].Last(s[0]),
          a1 = p.aut[1].Last(s[1])
        DO
          IF a0.symbol > a1.symbol THEN
            RETURN Arc{symbol := a0.symbol, dest := State{a0.dest, RedNull}}
          ELSIF a0.symbol < a1.symbol THEN
            RETURN Arc{symbol := a1.symbol, dest := State{RedNull, a1.dest}}
          ELSE
            RETURN Arc{symbol := a0.symbol, dest := State{a0.dest, a1.dest}}
          END
        END
      END;
    END;
  END Last;

PROCEDURE Rest(p: T; s: State): State =
  BEGIN
    IF s[0] = RedNull OR s[0] = RedUnit THEN
      IF s[1] = RedNull OR s[1] = RedUnit THEN
        <* ASSERT FALSE *>
      ELSE
        RETURN State{s[0], p.aut[1].Rest(s[1])}
      END;
    ELSE
      IF s[1] = RedNull OR s[1] = RedUnit THEN
        RETURN State{p.aut[0].Rest(s[0]), s[1]}
      ELSE
        WITH
          a0 = p.aut[0].Last(s[0]),
          a1 = p.aut[1].Last(s[1])
        DO
          IF a0.symbol > a1.symbol THEN
            RETURN State{p.aut[0].Rest(s[0]), s[1]}
          ELSIF a0.symbol < a1.symbol THEN
            RETURN State{s[0], p.aut[1].Rest(s[1])}
          ELSE
            RETURN State{p.aut[0].Rest(s[0]), p.aut[1].Rest(s[1])}
          END
        END
      END;
    END;
  END Rest;

PROCEDURE Root(p: T): State =
  BEGIN
    RETURN State{p.aut[0].Root(), p.aut[1].Root()}
  END Root;

(* DERIVED METHODS *)

PROCEDURE Final(p: T; s: State): Bools =
  BEGIN
    RETURN Bools{p.aut[0].Final(s[0]), p.aut[1].Final(s[1])}
  END Final;

PROCEDURE OutDeg(p: T; s: State): NAT =
  VAR deg: NAT := 0;
  BEGIN
    WHILE p.HasArcs(s) DO
      INC(deg);
      s := p.Rest(s)
    END;
    RETURN deg
  END OutDeg;

PROCEDURE First(p: T; s: State): Arc =
  VAR t: State := NullState;
  BEGIN
    WHILE p.HasArcs(s) DO
      t := s;
      s := p.Rest(s)
    END;
    IF t = NullState THEN <* ASSERT FALSE *> END;
    RETURN p.Last(t)
  END First;

PROCEDURE R(p: T; s: State; symbol: Symbol): State =
  BEGIN
    RETURN State{p.aut[0].R(s[0], symbol), p.aut[1].R(s[1], symbol)}
  END R;

PROCEDURE Walk(p: T; s: State; READONLY w: String): State =
  BEGIN
    RETURN State{p.aut[0].Walk(s[0], w), p.aut[1].Walk(s[1], w)}
  END Walk;

PROCEDURE EnumArcs(p: T; s: State; action: ArcAction) RAISES {Abort} =
  VAR i: NAT := 0;

  PROCEDURE DoEnum(r: State) RAISES {Skip, Abort} =
  (*
    Does EnumArcs on a given prefix "r" of the arcs out of "s".
    Raises Skip, Abort iff "action" does so. *)
    BEGIN
      IF NOT p.HasArcs(r) THEN
        (* Ok *)
      ELSE
        DoEnum(p.Rest(r));
        WITH a = p.Last(r) DO
          action(i, a);
          INC(i)
        END
      END
    END DoEnum;

  BEGIN
    TRY DoEnum(s) EXCEPT Skip => (* Ok *) END;
  END EnumArcs;

PROCEDURE EnumPaths(
    p: T;
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

    VAR final: Bools; (* p.Final(t); set by "EnumRest" at the end of the recursion. *)
        i: NAT := 0; (* Arc index (not including the NullLetter, if any) *)

    PROCEDURE EnumRest(r: State) RAISES {Skip, Abort} =
    (*
      Calls "enter" on "t" and enumerates a given prefix "r" of the arcs out of
   "t".
      Also sets the local variable "final" of DoEnumStates(t)) to p.Final(t).
      Raises Skip iff "enter" or "pop" raised "Skip". *)

      BEGIN
        IF NOT p.HasArcs(r) THEN
          final := p.Final(r);
          IF enter # NIL THEN enter(len, t, final) END;
        ELSE
          EnumRest(p.Rest(r));
          WITH a = p.Last(r) DO
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
    p: T;
    s: State;
    enter: StringAction := NIL;
    exit: StringAction := NIL;
  ) RAISES {Abort} =

  VAR rw: REF String := NEW(REF String, 100);

  PROCEDURE PathEnter (* : StateAction *) (
      len: NAT;
      s: State;
      final: Bools;
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
      final: Bools;
    ) RAISES {Skip, Abort} =
    BEGIN
      IF exit # NIL THEN exit(SUBARRAY(rw^, 0, len), s, final) END
    END PathExit;

  BEGIN
    p.EnumPaths(
      s,
      enter := PathEnter,
      push := PathPush,
      pop := PathPop,
      exit := PathExit
    );
  END EnumStrings;

PROCEDURE Accepts(p: T; s: State; op: BoolOp; READONLY w: String): BOOL =
  BEGIN
    CASE op OF
    | BoolOp.Union => RETURN p.aut[0].Accepts(s[0], w) OR p.aut[1].Accepts(s[1],
    w);
    | BoolOp.Inter => RETURN p.aut[0].Accepts(s[0], w) AND p.aut[1].Accepts(s[1], w);
    | BoolOp.Diff  => RETURN p.aut[0].Accepts(s[0], w) AND NOT p.aut[1].Accepts(s[1], w)
    | BoolOp.Symm  => RETURN p.aut[0].Accepts(s[0], w) # p.aut[1].Accepts(s[1], w);
    END;
  END Accepts;

PROCEDURE EnumSuffs(p: T; s: State; op: BoolOp; action: SuffixAction) RAISES {Abort} =

  PROCEDURE SuffEnter (* : StringAction *) (
      READONLY w: String;
      <*UNUSED*> d: State;
      final: Bools;
    ) RAISES {Skip,Abort} =
    VAR fb: BOOL;
    BEGIN
      (* Decide if string is in Suff(s[0]) op Suff(s[1]): *)
      CASE op OF
      | BoolOp.Union =>  fb := final[0] OR final[1];
      | BoolOp.Inter =>  fb := final[0] AND final[1];
      | BoolOp.Diff  =>  fb := final[0] AND NOT final[1];
      | BoolOp.Symm  =>  fb := final[0] # final[1];
      END;
      IF fb THEN action(w) END;
      (* Decide if it is worth continuing: *)
      CASE op  OF
      | BoolOp.Union => (* Ok *)
      | BoolOp.Inter => IF s[0] = RedNull OR s[1] = RedNull THEN RAISE Skip END;
      | BoolOp.Diff  => IF s[0] = RedNull THEN RAISE Skip END;
      | BoolOp.Symm  => (* Ok *)
      END
    END SuffEnter;

  BEGIN
    p.EnumStrings(s, enter := SuffEnter);
  END EnumSuffs;

PROCEDURE NSuffs(p: T; s: State; op: BoolOp; limit: NAT := LAST(NAT)): NAT RAISES {Abort} =

  VAR count: NAT :=0;

  PROCEDURE IncrCount(n: NAT) RAISES {Abort} =
    BEGIN
      IF count >= limit - n THEN RAISE Abort END;
      INC(count, n);
    END IncrCount;

  PROCEDURE TestCount(est: NAT) RAISES {Abort} =
  (*
    Given a lower bound "est" to the number of strings remaining,
    raises "Abort" if the "limit" is going to be reached or exceeded. *)
    BEGIN
      IF count >= limit - est THEN RAISE Abort END;
    END TestCount;

  PROCEDURE CountEnter (* : StateAction *) (
      <*UNUSED*> len: NAT;
      s: State;
      final: Bools;
    )  RAISES {Skip,Abort} =
    BEGIN
      (* Shortcuts: *)
      WITH a0 = p.aut[0], a1 = p.aut[1] DO
        IF s[0] = RedNull THEN
          CASE op OF
          | BoolOp.Union => IncrCount(a1.NSuffs(s[1])); RAISE Skip
          | BoolOp.Inter => RAISE Skip
          | BoolOp.Diff  => RAISE Skip
          | BoolOp.Symm  => IncrCount(a1.NSuffs(s[1])); RAISE Skip
          END;
        ELSIF s[1] = RedNull THEN
          CASE op OF
          | BoolOp.Union => IncrCount(a0.NSuffs(s[0])); RAISE Skip
          | BoolOp.Inter => RAISE Skip
          | BoolOp.Diff  => IncrCount(a0.NSuffs(s[0])); RAISE Skip
          | BoolOp.Symm  => IncrCount(a0.NSuffs(s[0])); RAISE Skip
          END;
        ELSE
          (* Do some quick estimates to see if we are going to exceed the limit: *)
          CASE op OF
          | BoolOp.Union => TestCount(MAX(a0.NSuffs(s[0]), a1.NSuffs(s[1])))
          | BoolOp.Inter => (* Ok *)
          | BoolOp.Diff  => TestCount(MAX(0, a0.NSuffs(s[0]) - a1.NSuffs(s[1])))
          | BoolOp.Symm  => TestCount(ABS(a0.NSuffs(s[0]) - a1.NSuffs(s[1])))
          END;
          (* Account for the current string, if in the counted set: *)
          WITH fm0 = final[0], fm1 = final[1]  DO
            CASE op  OF
            | BoolOp.Union =>  IF fm0 OR fm1 THEN IncrCount(1) END
            | BoolOp.Inter =>  IF fm0 AND fm1 THEN IncrCount(1) END
            | BoolOp.Diff  =>  IF fm0 AND NOT fm1 THEN IncrCount(1) END
            | BoolOp.Symm  =>  IF fm0 # fm1 THEN IncrCount(1) END
            END
          END;
          (* Continue the enumeration *)
        END;
      END;
    END CountEnter;

  BEGIN
    p.EnumPaths(s, enter := CountEnter);
    RETURN count
  END NSuffs;

PROCEDURE PrintSuffs(
    p: T;
    s: State;
    op: BoolOp;
    spr: StringPrinter.T;
  ) =

  PROCEDURE PrintSuffix (* : SuffixAction *) (READONLY w: String) RAISES {Abort} =
    BEGIN
      spr.PutString(w, rev := FALSE)
    END PrintSuffix;

  BEGIN
    TRY p.EnumSuffs(s, op, action := PrintSuffix) EXCEPT Abort => (* Ok *) END;
  END PrintSuffs;

BEGIN
END ReducedPair.

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
