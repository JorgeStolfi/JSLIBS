GENERIC INTERFACE GenSort(X, Y);

(* Generic sorting routine *)
(* See the copyright and disclaimer note at the end of this file. *)

PROCEDURE Sort(
    VAR (*IO*) s: ARRAY OF X.T; 
    VAR (*IO*) ss: REF ARRAY OF Y.T;
  );
  (* 
    QuickSort sorts  the elements in  the array "s".  If  "ss#NIL" then
    the elements  of "ss"  are exchanged in  parallel  with those  in "s".
    Formal parameter interface X must provide at least the type "X.T" and
    the comparison procedure "X.Compare" with the signature:

         PROCEDURE(a,b: X.T): [-1..+1]

    The elements of "s" will be sorted accordingly to the results of "X.Compare". *)

END GenSort.

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
