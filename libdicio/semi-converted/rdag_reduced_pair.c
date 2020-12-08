#ifndef _H
#define _H


/* See the copyright and disclaimer note at the end of this file. */

IMPORT Basics, StringPrinter, Reduced;
FROM Basics IMPORT NAT, BOOL, Skip, Abort;
FROM Reduced IMPORT Symbol, String;

REVEAL
  T == Public BRANDED OBJECT
    OVERRIDES

      HasArcs = HasArcs;
      Last = Last;
      Rest = Rest;
      Root = Root;
      Final = Final;

      OutDeg = OutDeg;
      First = First;
      R = R;
      Walk = Walk;

      EnumArcs = EnumArcs;
      EnumPaths = EnumPaths;
      EnumStrings = EnumStrings;

      Accepts = Accepts;
      EnumSuffs = EnumSuffs;
      NSuffs = NSuffs;
      PrintSuffs = PrintSuffs;

    ;};

CONST
  RedNull == Reduced.NullState;
  RedUnit == Reduced.UnitState;

PROCEDURE New(a0, a1: Reduced.T): T ==
  {
    return NEW(T, aut = ARRAY Which OF Reduced.T{a0, a1})
  ;} New;

/* PRIMITIVE METHODS */

PROCEDURE HasArcs(<*UNUSED); p: T; s: State): BOOL ==
  {
    return (s[0]!=RedNull)  AND  AND  (s[0]!=RedUnit)
       ) || ((s[1]!=RedNull)  AND  AND  (s[1]!=RedUnit)
  ;} HasArcs;

PROCEDURE Last(p: T; s: State): Arc ==
  {
    if ((s[0] == RedNull) || (s[0] == RedUnit)){
      if ((s[1] == RedNull) || (s[1] == RedUnit)){
        assert(FALSE );
      }else{
        with (a1 == p.aut[1].Last(s[1])){
          return Arc{symbol = a1.symbol, dest = State{RedNull, a1.dest}}
        ;}
      ;};
    }else{
      if ((s[1] == RedNull) || (s[1] == RedUnit)){
        with (a0 == p.aut[0].Last(s[0])){
          return Arc{symbol = a0.symbol, dest = State{a0.dest, RedNull}}
        ;}
      }else{
        with (
          a0 == p.aut[0].Last(s[0]),
          a1 == p.aut[1].Last(s[1])
       ){
          if ((a0.symbol > a1.symbol)){
            return Arc{symbol = a0.symbol, dest = State{a0.dest, RedNull}}
          }else if ((a0.symbol < a1.symbol)){
            return Arc{symbol = a1.symbol, dest = State{RedNull, a1.dest}}
          }else{
            return Arc{symbol = a0.symbol, dest = State{a0.dest, a1.dest}}
          ;}
        ;}
      ;};
    ;};
  ;} Last;

PROCEDURE Rest(p: T; s: State): State ==
  {
    if ((s[0] == RedNull) || (s[0] == RedUnit)){
      if ((s[1] == RedNull) || (s[1] == RedUnit)){
        assert(FALSE );
      }else{
        return State{s[0], p.aut[1].Rest(s[1])}
      ;};
    }else{
      if ((s[1] == RedNull) || (s[1] == RedUnit)){
        return State{p.aut[0].Rest(s[0]), s[1]}
      }else{
        with (
          a0 == p.aut[0].Last(s[0]),
          a1 == p.aut[1].Last(s[1])
       ){
          if ((a0.symbol > a1.symbol)){
            return State{p.aut[0].Rest(s[0]), s[1]}
          }else if ((a0.symbol < a1.symbol)){
            return State{s[0], p.aut[1].Rest(s[1])}
          }else{
            return State{p.aut[0].Rest(s[0]), p.aut[1].Rest(s[1])}
          ;}
        ;}
      ;};
    ;};
  ;} Rest;

PROCEDURE Root(p: T): State ==
  {
    return State{p.aut[0].Root(), p.aut[1].Root()}
  ;} Root;

/* DERIVED METHODS */

PROCEDURE Final(p: T; s: State): Bools ==
  {
    return Bools{p.aut[0].Final(s[0]), p.aut[1].Final(s[1])}
  ;} Final;

PROCEDURE OutDeg(p: T; s: State): NAT ==
  VAR deg: NAT = 0;
  {
    while (p.HasArcs(s)){
      INC(deg);
      s = p.Rest(s)
    ;};
    return deg
  ;} OutDeg;

PROCEDURE First(p: T; s: State): Arc ==
  VAR t: State = NullState;
  {
    while (p.HasArcs(s)){
      t = s;
      s = p.Rest(s)
    ;};
    if ((t == NullState)){ assert(FALSE ); ;};
    return p.Last(t)
  ;} First;

PROCEDURE R(p: T; s: State; symbol: Symbol): State ==
  {
    return State{p.aut[0].R(s[0], symbol), p.aut[1].R(s[1], symbol)}
  ;} R;

PROCEDURE Walk(p: T; s: State; READONLY w: String): State ==
  {
    return State{p.aut[0].Walk(s[0], w), p.aut[1].Walk(s[1], w)}
  ;} Walk;

PROCEDURE EnumArcs(p: T; s: State; action: ArcAction) RAISES {Abort} ==
  VAR i: NAT = 0;

  PROCEDURE DoEnum(r: State) RAISES {Skip, Abort} ==
  /*
    Does EnumArcs on a given prefix "r" of the arcs out of "s".
    Raises Skip, Abort iff "action" does so. */
    {
      if ((NOT p.HasArcs(r))){
        /* Ok */
      }else{
        DoEnum(p.Rest(r));
        with (a == p.Last(r)){
          action(i, a);
          INC(i)
        ;}
      ;}
    ;} DoEnum;

  {
    TRY DoEnum(s) EXCEPT Skip ==> /* Ok */ ;};
  ;} EnumArcs;

PROCEDURE EnumPaths(
    T p;
    State s;
    enter: StateAction = NULL;
    push: PathAction = NULL;
    pop: PathAction = NULL;
    exit: StateAction = NULL;
  ) RAISES {Abort} ==

  VAR len: NAT = 0;

  PROCEDURE DoEnumPaths(t: State) RAISES {Abort} ==
  /*
    Does EnumPaths starting at "t", a generic sucessor
    of "s". Assumes "t" is not NullState. */

    Bools *final; /* p.Final(t); set by "EnumRest" at the end of the recursion. */
        i: NAT = 0; /* Arc index (not including the NullLetter, if any) */

    PROCEDURE EnumRest(r: State) RAISES {Skip, Abort} ==
    /*
      Calls "enter" on "t" and enumerates a given prefix "r" of the arcs out of
   "t".
      Also sets the local variable "final" of DoEnumStates(t)) to p.Final(t).
      Raises Skip iff "enter" or "pop" raised "Skip". */

      {
        if ((NOT p.HasArcs(r))){
          final = p.Final(r);
          if ((enter!=NULL)){ enter(len, t, final) ;};
        }else{
          EnumRest(p.Rest(r));
          with (a == p.Last(r)){
            TRY
              if ((push!=NULL)){ push(len, t, i, a) ;};
              INC(len);
              DoEnumPaths(a.dest);
              DEC(len);
            EXCEPT
              Skip ==> /* Ok */
            ;};
            if ((pop!=NULL)){ pop(len, t, i, a) ;};
            INC(i);
          ;}
        ;}
      ;} EnumRest;

    {
      assert(t!=NullState );
      TRY EnumRest(t) EXCEPT Skip ==> /* Ok */ ;};
      if ((exit!=NULL)){
        <* FATAL Skip );
        {
          exit(len, t, final) /* Shouldn't raise "Skip" */
        ;}
      ;}
    ;} DoEnumPaths;

  {
    if ((s!=NullState)){ DoEnumPaths(s) ;};
  ;} EnumPaths;

PROCEDURE EnumStrings(
    T p;
    State s;
    enter: StringAction = NULL;
    exit: StringAction = NULL;
  ) RAISES {Abort} ==

  VAR rw: REF String = NEW(REF String, 100);

  PROCEDURE PathEnter /* : StateAction */ (
      NAT len;
      State s;
      Bools final;
    ) RAISES {Skip, Abort} ==
    {
      if ((enter!=NULL)){ enter(SUBARRAY(rw^, 0, len), s, final) ;}
    ;} PathEnter;

  PROCEDURE PathPush /* : PathAction */ (
      NAT len;
      <*UNUSED); org: State;
      <*UNUSED); i: NAT;
      Arc arc;
    ) RAISES {} ==
    {
      Basics.ExpandString(rw, len + 1);
      rw[len] = arc.symbol;
    ;} PathPush;

  PROCEDURE PathPop /* : PathAction */ (
      NAT len;
      <*UNUSED); org: State;
      <*UNUSED); i: NAT;
      Arc arc;
    ) RAISES {} ==
    {
      assert(rw[len] == arc.symbol );
    ;} PathPop;

  PROCEDURE PathExit /* : StateAction */ (
      NAT len;
      State s;
      Bools final;
    ) RAISES {Skip, Abort} ==
    {
      if ((exit!=NULL)){ exit(SUBARRAY(rw^, 0, len), s, final) ;}
    ;} PathExit;

  {
    p.EnumPaths(
      s,
      enter = PathEnter,
      push = PathPush,
      pop = PathPop,
      exit = PathExit
    );
  ;} EnumStrings;

PROCEDURE Accepts(p: T; s: State; op: BoolOp; READONLY w: String): BOOL ==
  {
    CASE op OF
    | BoolOp.Union ==> return p.aut[0].Accepts(s[0], w)) || (p.aut[1].Accepts(s[1],
    w);
    | BoolOp.Inter ==> return p.aut[0].Accepts(s[0], w))  AND  AND  (p.aut[1].Accepts(s[1], w);
    | BoolOp.Diff  ==> return p.aut[0].Accepts(s[0], w))  AND  AND  (NOT p.aut[1].Accepts(s[1], w)
    | BoolOp.Symm  ==> return p.aut[0].Accepts(s[0], w)!=p.aut[1].Accepts(s[1], w);
    ;};
  ;} Accepts;

PROCEDURE EnumSuffs(p: T; s: State; op: BoolOp; action: SuffixAction) RAISES {Abort} ==

  PROCEDURE SuffEnter /* : StringAction */ (
      String READONLY w;
      <*UNUSED); d: State;
      Bools final;
    ) RAISES {Skip,Abort} ==
    BOOL *fb;
    {
      /* Decide if string is in Suff(s[0]) op Suff(s[1]): */
      CASE op OF
      | BoolOp.Union ==>  fb = final[0]) || (final[1];
      | BoolOp.Inter ==>  fb = final[0])  AND  AND  (final[1];
      | BoolOp.Diff  ==>  fb = final[0])  AND  AND  (NOT final[1];
      | BoolOp.Symm  ==>  fb = final[0]!=final[1];
      ;};
      if ((fb)){ action(w) ;};
      /* Decide if it is worth continuing: */
      CASE op  OF
      | BoolOp.Union ==> /* Ok */
      | BoolOp.Inter ==> if ((s[0] == RedNull) || (s[1] == RedNull)){ RAISE Skip ;};
      | BoolOp.Diff  ==> if ((s[0] == RedNull)){ RAISE Skip ;};
      | BoolOp.Symm  ==> /* Ok */
      ;}
    ;} SuffEnter;

  {
    p.EnumStrings(s, enter = SuffEnter);
  ;} EnumSuffs;

PROCEDURE NSuffs(p: T; s: State; op: BoolOp; limit: NAT = LAST(NAT)): NAT RAISES {Abort} ==

  VAR count: NAT =0;

  PROCEDURE IncrCount(n: NAT) RAISES {Abort} ==
    {
      if ((count >= limit - n)){ RAISE Abort ;};
      INC(count, n);
    ;} IncrCount;

  PROCEDURE TestCount(est: NAT) RAISES {Abort} ==
  /*
    Given a lower bound "est" to the number of strings remaining,
    raises "Abort" if the "limit" is going to be reached or exceeded. */
    {
      if ((count >= limit - est)){ RAISE Abort ;};
    ;} TestCount;

  PROCEDURE CountEnter /* : StateAction */ (
      <*UNUSED); len: NAT;
      State s;
      Bools final;
    )  RAISES {Skip,Abort} ==
    {
      /* Shortcuts: */
      with (a0 == p.aut[0], a1 == p.aut[1]){
        if ((s[0] == RedNull)){
          CASE op OF
          | BoolOp.Union ==> IncrCount(a1.NSuffs(s[1])); RAISE Skip
          | BoolOp.Inter ==> RAISE Skip
          | BoolOp.Diff  ==> RAISE Skip
          | BoolOp.Symm  ==> IncrCount(a1.NSuffs(s[1])); RAISE Skip
          ;};
        }else if ((s[1] == RedNull)){
          CASE op OF
          | BoolOp.Union ==> IncrCount(a0.NSuffs(s[0])); RAISE Skip
          | BoolOp.Inter ==> RAISE Skip
          | BoolOp.Diff  ==> IncrCount(a0.NSuffs(s[0])); RAISE Skip
          | BoolOp.Symm  ==> IncrCount(a0.NSuffs(s[0])); RAISE Skip
          ;};
        }else{
          /* Do some quick estimates to see if we are going to exceed the limit: */
          CASE op OF
          | BoolOp.Union ==> TestCount(MAX(a0.NSuffs(s[0]), a1.NSuffs(s[1])))
          | BoolOp.Inter ==> /* Ok */
          | BoolOp.Diff  ==> TestCount(MAX(0, a0.NSuffs(s[0]) - a1.NSuffs(s[1])))
          | BoolOp.Symm  ==> TestCount(ABS(a0.NSuffs(s[0]) - a1.NSuffs(s[1])))
          ;};
          /* Account for the current string, if in the counted set: */
          with (fm0 == final[0], fm1 == final[1] ){
            CASE op  OF
            | BoolOp.Union ==>  if ((fm0) || (fm1)){ IncrCount(1) ;}
            | BoolOp.Inter ==>  if ((fm0)  AND  AND  (fm1)){ IncrCount(1) ;}
            | BoolOp.Diff  ==>  if ((fm0)  AND  AND  (NOT fm1)){ IncrCount(1) ;}
            | BoolOp.Symm  ==>  if ((fm0!=fm1)){ IncrCount(1) ;}
            ;}
          ;};
          /* Continue the enumeration */
        ;};
      ;};
    ;} CountEnter;

  {
    p.EnumPaths(s, enter = CountEnter);
    return count
  ;} NSuffs;

PROCEDURE PrintSuffs(
    T p;
    State s;
    BoolOp op;
    spr: StringPrinter.T;
  ) ==

  PROCEDURE PrintSuffix /* : SuffixAction */ (READONLY w: String) RAISES {Abort} ==
    {
      spr.PutString(w, rev = FALSE)
    ;} PrintSuffix;

  {
    TRY p.EnumSuffs(s, op, action = PrintSuffix) EXCEPT Abort ==> /* Ok */ ;};
  ;} PrintSuffs;

{
;} ReducedPair.

/****************************************************************************/
/* (C) Copyright 1992 Universidade Estadual de Campinas (UNICAMP)           */
/*                    Campinas, SP, Brazil                                  */
/*                                                                          */
/* Authors:                                                                 */
/*                                                                          */
/*   Tomasz Kowaltowski  - CS Dept, UNICAMP <tomasz@dcc.unicamp.br>         */
/*   Claudio L. Lucchesi - CS Dept, UNICAMP <lucchesi@dcc.unicamp.br>       */
/*   Jorge Stolfi        - CS Dept, UNICAMP <stolfi@dcc.unicamp.br>         */
/*                                                                          */
/* This file can be freely distributed, modified, and used for any          */
/*   non-commercial purpose, provided that this copyright and authorship    */
/*   notice be included in any copy or derived version of this file.        */
/*                                                                          */
/* DISCLAIMER: This software is offered ``as is'', without any guarantee    */
/*   as to fitness for any particular purpose.  Neither the copyright       */
/*   holder nor the authors or their employers can be held responsible for  */
/*   any damages that may result from its use.                              */
/****************************************************************************/
