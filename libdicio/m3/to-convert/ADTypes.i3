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

INTERFACE ADTypes;

TYPE
  INT = INTEGER;
  NAT = CARDINAL;
  BOOL = BOOLEAN;
  
  IntPoint = ARRAY [0..1] OF INT;     (* Indices of grid point *)

  RealPoint = ARRAY [0..1] OF REAL;   (* Cartesian coordinates *)

  RealRange = RECORD lo, hi: REAL END;  (* Set OF REALs between /lo/ and /hi/, incl. *)

  IntRange = RECORD lo, hi: INT END;    (* Set of INTs between /lo/ and /hi/, incl. *)
  
  IntBox = ARRAY [0..1] OF IntRange;    (* Horizontal and vertical ranges *)

  Dims = ARRAY [0..1] OF REAL;       (* Pair of dimensions *)

  Counts = ARRAY [0..1] OF NAT;      (* Pair of counts *)

TYPE 
  ArcControlPoints = ARRAY [0..4] OF RealPoint; 
    (* 
      The arc with control points "cp" consists of a straight segment
      from cp[0] to cp[1], followed by a Bezier cubic arc with control
      points cp[1..4]. *)

CONST
  Nowhere = IntPoint{LAST(NAT), LAST(NAT)};

PROCEDURE ClipIntPoint(gp: IntPoint; box: IntBox): IntPoint;
  (*
    Returns the point in "box" that is nearest to "gp". *)

END ADTypes.
