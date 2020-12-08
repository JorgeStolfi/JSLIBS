INTERFACE SortedReduced;

(* Reduced deterministic acyclic automata with sorted transitions. *)
(* See the copyright and disclaimer note at the end of this file. *)

(*
  A SortedReduced.T is a special case of Reduced.T such that 
  "Last(s)" is always the transition out ot "s" with largest symbol
  label.  
  
  This convention allows slightly faster rejection of words, since the
  "R" method can establish the absence of a transition without
  scanning the entire list.  It also allows much faster
  addition/deletion, since state uniqueness in the Reduced.T sense
  becomes a corollary of DAG state uniqueness --- which means
  the "SetArc" method does not need to do costly set hashing and 
  comparison.

  On the other hand, the ordering of transition prevents efficient
  sharing of substates; thus, the "Fold" method is a no-op for these
  automata.
*)  

TYPE
  T <: Public;
  Public = Reduced.T OBJECT
    METHODS

      (*
        The folowing methods from "Reduced.T" have their semantics
        redefined in a "SortedReduced.T":
        
          "Last":  returns always the transition with highest symbol label.
          
          "Rest(s)":  returns the given state minus 
          
          "Fold()":  a no-op.
          
        The following derived methods have their semantics affected
        by these changes:
        
          "First(s)":  returns the transition with lowest symbol label.
          
          "R(s, x)":  the cost is "O(k+1)", where "k" is the numebr of arcs
             out of "s" that have labels greater than "x".
                         
          "SetArc(s, x, d)":  likely to be faster.  Its cost may be
            proportional to the number of arcs out of "s" with labels
            greater than "symbol"; hence it is a good idea to build states
            by adding the arcs in increasing symbol order.
          
          "SetFinal(s, b)":  likely to be faster.
          
          "Walk(s, w)":  likely to be faster if "w" is not in the language.
          
          "Accepts(s, w)":  ditto.
          
          "Rank(s, w)":  about twice as fast.
          
          "AddString", "SubString", and relatives:  somewhat slower.
          
          "EnumOutArcs":  the arcs are visited in order of increasing
             symbol label.
             
          "EnumInArcs":  arcs with same origin are ordered by increasing
             symbol label.
             
          "EnumPaths":  paths are listed in increasing lexicographic order.
                       
      *)

      AppendArc(s: State; symbol: Symbol; dest: State): State RAISES {Full};
      (*
        Returns the (possibly new) state that has the same final bit
        and outgoing arcs as "s", plus one new arc labeled with
        "symbol" and leading to "dest".  Requires "symbol" to be
        strictly greater than the label of any proper arc out of "s".

        The "dest" state can be NullState, in which case the
        arc is not added, and the result is just "s" itself. *)

    END;

END SortedReduced.

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
