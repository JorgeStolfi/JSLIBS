#ifndef _H
#define _H


/* Cartesian product of reduced automata */
/* See the copyright and disclaimer note at the end of this file. */

/*
  A ReducedPair.T "p" is the Cartesian product of two reduced automata,
  p.aut[0] and p.aut[1].

  Each state "s" of "p" is a pair of Reduced.States: one state s[0]
  from p.aut[0] and one state s[1] from p.aut[1].  The state is proper
  iff at least one of the component states is proper.

  There is a transition in "p" from state "s" to state "t" by the symbol
  "x" if and only if there is a transition from s[0] to t[0]
  in p.aut[0], and one from s[1] to t[1] in p.aut[1], both labeled "x".
  The transition is proper iff at least one of the transitions is
  proper; that is, at least one of s[0] and s[1] is non-null,
  and at least one of t[0] and t[1] is non-null.

  In particular, if p.aut[0] has a proper arc from s[0] to t[0]
  labeled "x", but p.aut[1] has no arc labeled "x" out of s[1],
  then "p" has a proper arc labeled "x" from "s" to the pair
  (t[0], NullState).

  Note that the virtual state ReducedPair.NullState
  (the state pair whose two components are both Reduced.NullState)
  is its own successor under every symbol.

  Obviously, each State of the product has two "final" flags,
  one from each Reduced.State.  Thus, a ReducedPair.T has four
  "naked" states (with no proper transitions): (N,N), (N,U), (U,N),
  and (U,U), where N stands for Reduced.NullState,
  and U for Reduced.UnitState.
  */

IMPORT Reduced, StringPrinter;
FROM Basics IMPORT NAT, BOOL, Skip, Abort;
FROM Reduced IMPORT Symbol, String;

TYPE
  Which == [0..1];                /* Which component of the pair? */

  State == ARRAY Which OF Reduced.State;          /* A state of the product automaton */
  Arc == RECORD symbol: Symbol; dest: State ;};  /* A transition of the prod. aut. */
  Bools == ARRAY Which OF BOOL;                   /* The "final" marks of a state. */

TYPE
  BoolOp == {
    Union, /* Set union */
    Inter, /* Set intersection */
    Diff,  /* Set difference, A\B */
    Symm   /* Symmetric set difference (XOR, A\B U B\A) */
  };

CONST
  NullState == State{Reduced.NullState, Reduced.NullState};

TYPE
  T <: Public;

  Public == OBJECT
      aut: ARRAY Which OF Reduced.T;
    METHODS

      HasArcs(s: State): BOOL;
      /*
        TRUE iff the state "s" has at least one proper outgoing arc. */

      Final(s: State): Bools;
      /*
        The "final" bits of s[0] in p.aut[0] and s[1] in p.aut[1]. */

      Last(s: State): Arc;
      /*
        The proper arc out of "s" with highest symbol label.
        Requires HasArcs(s). */

      Rest(s: State): State;
      /*
        The state "s" minus "Last(s)", that is, a state that has
        the same proper outgoing arcs and final bits as "s",
        except for the arc with highest symbol label.
        Requires HasArcs(s). */

      Root(): State;
      /*
        Returns the current root state (the pair of roots of the
        factor automata). */

      /******************************************************************/
      /* DERIVED METHODS                                                */
      /******************************************************************/

      /*
        These methods could be ordinary procedures, since
        they can be defined in terms of the primitives above;
        they are declared here as methods to reduce client confusion. (?) */

      OutDeg(s: State): NAT;
      /*
        Number of proper arcs out of state "s" (0 iff s == NullState 
        or s == UnitState).
        Cost: O(OutDeg(s)) time, 0 space. */

      First(s: State): Arc;
      /*
        The first arc out of state "s". Requires "s" to be non-null.
        Cost: O(OutDeg(s)) time, 0 space. */

      R(s: State; symbol: Symbol): State;
      /*
        The successor of the given state through the arc labeled with
        the given symbol.  Returns NullState if "s" is NullState or has
        no outgoing edge labeled with that "symbol". */

      Walk(s: State; READONLY w: String): State;
      /*
        Spells "w" in the automaton starting from "s", that is,
        returns the state reached from "s" by the unique path
        whose edges are labeled with the symbols of "w". */

      /******************************************************************/
      /* TRAVERSAL                                                      */
      /******************************************************************/

      EnumArcs(s: State; action: ArcAction) RAISES {Skip, Abort};
      /*
        Just like Reduced.EnumArcs, for the product automaton. */

      EnumPaths(
          State s;
          enter: StateAction = NULL;
          push: PathAction = NULL;
          pop: PathAction = NULL;
          exit: StateAction = NULL;
        ) RAISES {Abort};
      /*
        Just like Reduced.EnumPaths, but for the product automaton.
        Enumerates only paths consisting entirely of proper states
        and arcs.
        */

      EnumStrings(
          State s;
          enter: StringAction = NULL;
          exit: StringAction = NULL;
        ) RAISES {Abort};
      /*
        Like EnumPaths, but passes the whole String spelled by the
        current path to the "enter" and "exit" actions.
        */

      /*****************************************************************/
      /* BOOLEAN OPERATIONS ON SUFFIX SETS                             */
      /*****************************************************************/

      Accepts(s: State; op: BoolOp; READONLY w: String): BOOL;
      /*
        TRUE iff "w" is in the language "Suff(s[0]) op Suff(s[1])".
        Equivalent to "f[0] op f[1]" where f == Final(Walk(s, w)). */

      EnumSuffs(s: State; op: BoolOp; action: SuffixAction) RAISES {Abort};
      /*
        Enumerates the specified Boolean combination of the
        Suff(s[0]) in pair.aut[0] and Suff(s[1]) in pair.aut[1],
        in lexicographic Symbol order, applying "action" to each suffix. */

      NSuffs(s: State; op: BoolOp; limit: NAT = LAST(NAT)): NAT RAISES {Abort};
      /*
        Counts the number of strings in the specified Boolean
        combination of Suff(s[0]) and Suff(s[1]).
        If the count is "limit" or more, raises "Abort". */

      PrintSuffs(s: State; op: BoolOp; spr: StringPrinter.T);
      /*
        Prints the the specified Boolean combination of the suffix sets
        Suff(s[0]) in pair.aut[0] and Suff(s[1]) in pair.aut[1],
        in lexicographic Symbol order, formatting as specified
        by "spr". */

  ;};

TYPE
  ArcAction == PROCEDURE(i: NAT; arc: Arc) RAISES {Skip,Abort};
  PathAction == PROCEDURE(len: NAT; org: State; i: NAT; a: Arc) RAISES {Skip,Abort};
  StateAction == PROCEDURE(len: NAT; s: State; final: Bools) RAISES {Skip,Abort};
  StringAction == PROCEDURE(READONLY w: String; d: State; final: Bools) RAISES {Skip,Abort};
  SuffixAction == PROCEDURE(READONLY w: String) RAISES {Skip,Abort};

/* CREATION */

PROCEDURE New(a0, a1: Reduced.T): T;
/*
  The product of the automata "a0" and "a1". */

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
