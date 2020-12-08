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

MODULE ADMath;

IMPORT Math;

PROCEDURE Tan(u: REAL): REAL =
  BEGIN
    RETURN FLOAT(Math.tan(FLOAT(u, LONGREAL)))
  END Tan;

PROCEDURE Cos(u: REAL): REAL =
  BEGIN
    RETURN FLOAT(Math.cos(FLOAT(u, LONGREAL)))
  END Cos;

PROCEDURE Sin(u: REAL): REAL =
  BEGIN
    RETURN FLOAT(Math.sin(FLOAT(u, LONGREAL)))
  END Sin;
  
PROCEDURE Dir(x, y: REAL): REAL =
  BEGIN
    RETURN FLOAT(Math.atan2(FLOAT(y, LONGREAL), FLOAT(x, LONGREAL)))
  END Dir;

PROCEDURE Sqrt(u: REAL): REAL =
  BEGIN
    RETURN FLOAT(Math.sqrt(FLOAT(u, LONGREAL)))
  END Sqrt;

PROCEDURE Sqr(u: REAL): REAL =
  BEGIN
    RETURN u * u
  END Sqr;

PROCEDURE L(x: REAL): LONGREAL =
  BEGIN
    RETURN FLOAT(x, LONGREAL) (* Modula-3 is a crock. *)
  END L;

<*UNUSED*>
PROCEDURE F(x: LONGREAL): REAL =
  BEGIN
    RETURN FLOAT(x, REAL) (* Modula-3 is a crock. *)
  END F;

BEGIN
END ADMath.
