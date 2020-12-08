INTERFACE Basics;

(* Basic types: symbols, strings, and things like that. *)
(* See the copyright and disclaimer note at the end of this file. *)

IMPORT Word;

CONST
  SymbolBits = 8;
  MaxSymbol = Word.Minus(Word.Shift(1, SymbolBits), 1);

TYPE
  Symbol = BITS SymbolBits FOR [0..MaxSymbol];
  String = ARRAY OF Symbol;

TYPE
  INT = INTEGER;
  NAT = CARDINAL;
  POS = [1..LAST(CARDINAL)];
  BOOL = BOOLEAN;
  CHARS = ARRAY OF CHAR;

PROCEDURE ExpandString(VAR s: REF String; n: CARDINAL);
PROCEDURE ExpandChars(VAR s: REF CHARS; n: CARDINAL);
(*
  If "s^" has "n" or more elements, does nothing.
  If "s^" has less than "n" elements, allocates a new REF String
  or REF CHARS "r" with at least "n" elements,
  copies "s^" into "r^[0..LAST(s^)]", and sets "s := r".

  Useful when building things one element at a time. *)

PROCEDURE CopyString(READONLY s: String): REF String;
(*
  Returns a copy of "s" in a newly allocated REF String. *)

EXCEPTION
  Done;  (* Means that we reached the end of some stream or sequence. *)
  Full;  (* Means that some preallocated storage was exhausted. *)
  Skip;  (* To enumeration procedures, means "skip this part" *)
  Abort; (* To enumeration procedures, means "stop the enumeration". *)

END Basics.

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
