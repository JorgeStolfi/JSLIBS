INTERFACE ParamUtil;

(* Utilities for command line parsing/printing. *)
(* See the copyright and disclaimer note at the end of this file. *)

(*
  This interfaces provides some customized versions of the "ParseParams"
  methods. 

  The argument parsing routines ("GetText", "GetRd", etc.) all print
  the parsed options to the "wr" writer (if it is not NIL), using
  "PrintBoolArg", "PrintTextArg", etc.
*)

IMPORT ParseParams, Wr, Rd;
FROM Basics IMPORT INT, BOOL;

PROCEDURE GetText(pp: ParseParams.T; wr: Wr.T; key: TEXT; default: TEXT): TEXT
  RAISES {ParseParams.Error};
  (*
    Roughly equivalent to "IF pp.keywordPresent(key) 
    THEN RETURN pp.getNext() ELSE RETURN default END" *)

PROCEDURE GetBool(pp: ParseParams.T; wr: Wr.T; key: TEXT): BOOL 
  RAISES {};
  (*
    Roughly equivalent to "pp.keywordPresent(key)". *)

PROCEDURE GetInt(
    pp: ParseParams.T; 
    wr: Wr.T; 
    key: TEXT; 
    default: INT;
    min: INT := FIRST(INT);
    max: INT := LAST(INT);
  ): INT 
  RAISES {ParseParams.Error};
  (*
    Roughly equivalent to "IF pp.keywordPresent(key) 
    THEN RETURN pp.getNextInt(min, max) 
    ELSE RETURN default END" *)

PROCEDURE CheckFileName(pp: ParseParams.T; name:TEXT): TEXT 
  RAISES {ParseParams.Error};
  (*
    Raises "ParseParams.Error" if "name" is empty, or starts with "-"
    but is not just "-". Otherwise returns "name" unchanged. *)

PROCEDURE GetFileName(pp: ParseParams.T; wr: Wr.T; key: TEXT): TEXT 
  RAISES {ParseParams.Error};
  (*
    Parses a filename argument, preceded by "key".  Returns the empty
    string if "key" is not present in the command line. 
    Checks the filename for validity.
    Roughly equivalent to "IF pp.keywordPresent(key) 
    THEN RETURN CheckFileName(pp.getNext()) ELSE RETURN "\"\"" END". *)

PROCEDURE GetLogWr(pp: ParseParams.T; wr: Wr.T): Wr.T
  RAISES {ParseParams.Error};
  (*
    Parses the optional argument "-log <filename>". 
    (For compatibility, "-mess <filename>" is also acceptable.)

    Returns "Stdio.stderr", if the "-log" key is not present, or the
    filename given is "-"; otherwise opens and returns a writer to the
    specified file.

    Raises "ParseParams.Error" if the given filename is empty, or starts
    with "-" but is not just "-", or if the file cannot be opened. *)

PROCEDURE GetRd(pp: ParseParams.T; wr: Wr.T; key: TEXT): Rd.T
  RAISES {ParseParams.Error};
  (*
    If the "key" argument was given, opens and returns a reader to the
    file whose name follows it (or "Stdio.stdin" if the filename
    given is "-").  If the "key" was not given, returns NIL.

    Raises "ParseParams.Error" if the given filename is empty, or starts
    with "-" but is not just "-", or if the file cannot be opened. *)

PROCEDURE GetWr(pp: ParseParams.T; wr: Wr.T; key: TEXT): Wr.T
  RAISES {ParseParams.Error};
  (*
    If the "key" argument was given, opens and returns a writer to the
    file whose name follows it (or "Stdio.stdout" if the filename
    given is "-").  If the "key" was not given, returns NIL.

    Raises "ParseParams.Error" if the given filename is empty, or starts
    with "-" but is not just "-", or if the file cannot be opened. *)

PROCEDURE PrintTextArg(wr: Wr.T; key: TEXT; arg: TEXT);
PROCEDURE PrintBoolArg(wr: Wr.T; key: TEXT);
PROCEDURE PrintIntArg(wr: Wr.T; key: TEXT; arg: INT);
  (*
    Prints "key" and "arg" to "wr", in the format 
    used by other procedures in this interface. *)

END ParamUtil.

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
 
