INTERFACE StringPrinter;

(* Auto-formatted string printout. *)
(* See the copyright and disclaimer note at the end of this file. *)

(*
  A StringPrinter.T is an object for printing strings to a
  given writer.

  The StringPrinter takes care of symbol-to-character conversion,
  string separators, indentation, line breaking, and output truncation,
  acording to parameters given at creation time. *)

IMPORT Wr, Encoding;
FROM Basics IMPORT String, Abort;
FROM Basics IMPORT NAT, BOOL;

TYPE

  T <: Public;

  Public = OBJECT METHODS

      PutString(READONLY w: String; rev: BOOL := FALSE) RAISES {Abort};
      (*
        Prints "w" to the underlying writer, with the appropriate
        conversions, separators, and line breaks, as specified
        by the parameters given to "Init".

        If "rev" is TRUE, prints the symbols of "w" in reverse order.

        "PutString" will use the "encoding"'s "PrintLetter" method
        to output each Symbol of "w"; unless "w" has zero length, 
        in which case "PutString" will write the specified
        "empty" text.

        If "w" is not the first string printed since the
        last "Init" call, then, before writing "w",
        "PutString" will output the "sep" text.

        After writing the eventual separator, if writing "w" 
        would cause the total number of characters  printed to 
        exceed "maxChars", then "PutString" will output the 
        "etc" text instead of "w", and raise "Abort".

        Otherwise, if writing "w" would cause the current column 
        to exceed "rightMargin", "PutString" will output a newline character
        and "leftMargin" blanks.  For this purpose, "PutString" assumes
        that right after the "Init" call the writer is positioned at
        the "initialColumn" (counting from 0).

        The "maxChars" and "rightMargin" limits assume that one
        Basics.Symbol is printed by "lpr" as one character; they include the
        lengths of separators and empty strings, but do not
        include the "etc" text that gets printed when "maxChars"
        is exceeded.  Note also that, due to the order of tests,
        the "rightmargin" may be exceeded by Length(sep) + Length(etc).
        *)

      Reset();
      (*
        Re-initializes "self" with the same parameters 
        used in the last "Init" call.  In particular,
        resets the total character count to zero, and the 
        current column to the "initialColumn" given to "Init".
        *)

    END;

PROCEDURE Init(
    e: T;
    wr: Wr.T;
    encoding: Encoding.T;
    empty: TEXT := "()";
    sep: TEXT := ", ";
    maxChars: NAT := LAST(NAT);
    etc: TEXT := "...";
    rightMargin: NAT := 60;
    leftMargin: NAT := 0;
    initialColumn: NAT := 0;
  );
  (*
    Initializes "self" with the given parameters.
    For the meaning of the parameters, see the description of 
    the "PutString" method above.
    *)

PROCEDURE New(
    wr: Wr.T;
    encoding: Encoding.T;
    empty: TEXT := "()";
    sep: TEXT := ", ";
    maxChars: NAT := LAST(NAT);
    etc: TEXT := "...";
    rightMargin: NAT := 60;
    leftMargin: NAT := 0;
    initialColumn: NAT := 0;
  ): T;

END StringPrinter.

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
