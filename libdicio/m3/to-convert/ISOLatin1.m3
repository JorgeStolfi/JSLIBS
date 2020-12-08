(* Last edited on 1999-06-05 20:32:58 by stolfi *)

MODULE ISOLatin1;

(*
  Based on "ISOChar.m3" by Jim Meehan and Henri Gouraud 
  (Copyright 1994, Digital Equipment Corp.)
*)

IMPORT ASCII, Word;

BEGIN
  FOR c := '\000' TO '\377' DO Upper [c] := ASCII.Upper [c] END;
  FOR c := '\340' TO '\376' DO
    IF c # '\367' THEN
      Upper [c] := VAL (ORD (c) - ORD ('a') + ORD ('A'), CHAR)
    END
  END;

  FOR c := '\000' TO '\377' DO Lower [c] := ASCII.Lower [c] END;
  FOR c := '\300' TO '\336' DO
    IF c # '\327' THEN
      Lower [c] := VAL (ORD (c) - ORD ('A') + ORD ('a'), CHAR)
    END
  END;

  FOR c := '\000' TO '\377' DO
    IF c IN Graphics THEN
      Control [c] := VAL (Word.And (ORD (c), 8_37), CHAR)
    ELSE
      Control [c] := c
    END
  END
END ISOLatin1.

(****************************************************************************)
(* (C) Copyright 1995 Universidade Estadual de Campinas (UNICAMP)           *)
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
