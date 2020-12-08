(****************************************************************************)
(* (C) Copyright 1992 Universidade Estadual de Campinas (UNICAMP)           *)
(*                    Campinas, SP, Brazil                                  *)
(*                                                                          *)
(* Authors:                                                                 *)
(*                                                                          *)
(*   Tomasz Kowaltowski  - CS Dept, UNICAMP <tomasz@dcc.unicamp.ansp.br>    *)
(*   Claudio L. Lucchesi - CS Dept, UNICAMP <lucchesi@dcc.unicamp.ansp.br>  *)
(*   Jorge Stolfi        - DEC Systems Research Center <stolfi@src.dec.com> *)
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

INTERFACE Adviser;
  
(* HISTORY:
   26/Aug/92: Adapted for Modula-3 version-2.04; minor modifications
*)   
  
(* Advising module for the Portuguese language, based on its nearly 
phonetic spelling.
*)
  
  IMPORT ReadOnlyPermDAG, AVLTextTree;
  FROM Basics IMPORT Abort;

  PROCEDURE Alternatives(READONLY w: TEXT;
                         Aut: ReadOnlyPermDAG.T;
                         LocalTree: AVLTextTree.T;
                         VAR (*IO*) AltTree: AVLTextTree.T;
                         Transpose,
                         Remove,
                         Subst,
                         Insert,
                         ForceAlternatives: BOOLEAN) RAISES{Abort};
                         
  (* Given a word /w/, an automaton /aut/,  and  an AVLTree /LocalTree/,
  the procedure returns in  the tree /AltTree/ alternative spellings for
  /w/.   /Transpose/, /Remove/, /Subst/  and /Insert/ force transposing,
  removal, substitution and  insertion of letters,  for words of minimum
  lengths.  If  /ForceAlternative/ is true  then  all  methods are tried
  even for short words. *)
 

END Adviser.
