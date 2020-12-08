(****************************************************************************)
(* (C) Copyright 1992 Universidade Estadual de Campinas (UNICAMP)           *)
(*                    Campinas, SP, Brazil                                  *)
(*                                                                          *)
(* Authors:                                                                 *)
(*                                                                          *)
(*   Tomasz Kowaltowski  - CS Dept, UNICAMP <tomasz@dcc.unicamp.ansp.br>    *)
(*   Claudio L. Lucchesi - CS Dept, UNICAMP <lucchesi@dcc.unicamp.ansp.br>  *)
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

INTERFACE ADMath;

PROCEDURE Tan(u: REAL): REAL;

PROCEDURE Cos(u: REAL): REAL;
  
PROCEDURE Sin(u: REAL): REAL;
  
PROCEDURE Dir(x, y: REAL): REAL;

PROCEDURE Sqrt(u: REAL): REAL;

PROCEDURE Sqr(u: REAL): REAL;    (* = u*u *)

PROCEDURE L(x: REAL): LONGREAL;  (* = FLOAT(x, LONGREAL) *)

PROCEDURE F(x: LONGREAL): REAL;  (* = FLOAT(x, REAL) *)

END ADMath.


