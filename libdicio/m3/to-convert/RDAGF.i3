INTERFACE RDAGF;

(* Implementation details and low-level operations for "RDAG.T" *)
(* See the copyright and disclaimer note at the end of this file. *)

TYPE
  T = RDAG.T;

REVEAL
  RDAG.Parent = DAG.T BRANDED OBJECT END;
    (*
      An "RDAGF.T" is implemented as a "DAG.T" with a linked hash table 
      that allows one to find quickly a node given its data fields.

      The hash table is rebuilt (lazily) whenever there are changes in the
      epoch or the allocated size.  It needs to be completed whenever
      new nodes are added.

    *)

      
  




END RDAGF.

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
