MODULE ParamUtil;

(* See the copyright and disclaimer note at the end of this file. *)

IMPORT ParseParams, Wr, Rd, Text, Fmt, Thread, Util;
FROM Basics IMPORT INT, BOOL;
FROM Stdio IMPORT stderr;

PROCEDURE KeywordPresent(pp: ParseParams.T; key: TEXT): BOOL =
  BEGIN
    RETURN pp.keywordPresent(key) OR pp.keywordPresent(Util.ToLowerCase(key))
  END KeywordPresent;

PROCEDURE GetText(pp: ParseParams.T; wr: Wr.T; key: TEXT; default: TEXT): TEXT
  RAISES {ParseParams.Error} =
  BEGIN
    IF KeywordPresent(pp, key) THEN
      WITH r = pp.getNext() DO
        IF wr # NIL THEN PrintTextArg(wr, key, r) END;  
        RETURN r
      END
    ELSE
      RETURN default
    END
  END GetText;

PROCEDURE GetBool(pp: ParseParams.T; wr: Wr.T; key: TEXT): BOOL 
  RAISES {} =
  BEGIN
    IF KeywordPresent(pp, key) THEN
      IF wr # NIL THEN PrintBoolArg(wr, key) END;  
      RETURN TRUE
    ELSE
      RETURN FALSE
    END
  END GetBool;

PROCEDURE GetInt(
    pp: ParseParams.T; 
    wr: Wr.T; 
    key: TEXT; 
    default: INT;
    min: INT := FIRST(INT);
    max: INT := LAST(INT);
  ): INT 
  RAISES {ParseParams.Error} =
  BEGIN
    IF KeywordPresent(pp, key) THEN
      WITH r = pp.getNextInt(min, max) DO
        IF wr # NIL THEN PrintIntArg(wr, key, r) END;  
        RETURN r
      END
    ELSE
      RETURN default
    END
  END GetInt;

PROCEDURE GetFileName(pp: ParseParams.T; wr: Wr.T; key: TEXT): TEXT 
  RAISES {ParseParams.Error} =
  BEGIN
    IF KeywordPresent(pp, key) THEN
      WITH fileName = pp.getNext() DO
        IF wr # NIL THEN PrintTextArg(wr, key, fileName) END;
        RETURN CheckFileName(pp, fileName)
      END
    ELSE
      RETURN ""
    END;
  END GetFileName;

PROCEDURE CheckFileName(pp: ParseParams.T; name: TEXT): TEXT 
  RAISES {ParseParams.Error} =
  BEGIN
    WITH len = Text.Length(name) DO
      IF len = 0 THEN 
        pp.error("file name should not be empty")
      ELSIF len > 1 AND Text.GetChar(name, 0) = '-' THEN 
        pp.error("file name should not begin with \"-\"")
      END;
      RETURN name
    END;
  END CheckFileName;

PROCEDURE GetLogWr(pp: ParseParams.T; wr: Wr.T): Wr.T
  RAISES {ParseParams.Error} =
  BEGIN
    IF pp.keywordPresent("-log") 
    OR pp.keywordPresent("-mess") 
    THEN
      WITH fileName = pp.getNext() DO
        IF wr # NIL THEN PrintTextArg(wr, "-log", fileName) END;
        EVAL CheckFileName(pp, fileName);
        IF Text.Equal(fileName, "-") THEN
          RETURN stderr
        ELSE
          RETURN Util.OpenWr(fileName, wr)
        END
      END
    ELSE
      RETURN stderr
    END
  END GetLogWr;

PROCEDURE GetRd(pp: ParseParams.T; wr: Wr.T; key: TEXT): Rd.T
  RAISES {ParseParams.Error} =
  BEGIN
    WITH fileName = GetFileName(pp, wr, key) DO
      RETURN Util.OpenRd(fileName, wr)
    END
  END GetRd;

PROCEDURE GetWr(pp: ParseParams.T; wr: Wr.T; key: TEXT): Wr.T
  RAISES {ParseParams.Error} =
  BEGIN
    WITH fileName = GetFileName(pp, wr, key) DO
      RETURN Util.OpenWr(fileName, wr)
    END
  END GetWr;

PROCEDURE PrintTextArg(wr: Wr.T; key: TEXT; arg: TEXT) =
  <* FATAL Thread.Alerted, Wr.Failure *>
  BEGIN
    Wr.PutText(wr, " \\\n");
    Wr.PutText(wr, "  ");
    Wr.PutText(wr, key);
    Wr.PutText(wr, " ");
    DoPrintTextArg(wr, arg);
  END PrintTextArg;

PROCEDURE PrintBoolArg(wr: Wr.T; key: TEXT) =
  <* FATAL Thread.Alerted, Wr.Failure *>
  BEGIN
    Wr.PutText(wr, " \\\n");
    Wr.PutText(wr, "  ");
    Wr.PutText(wr, key);
  END PrintBoolArg;

PROCEDURE PrintIntArg(wr: Wr.T; key: TEXT; arg: INT) =
  <* FATAL Thread.Alerted, Wr.Failure *>
  BEGIN
    Wr.PutText(wr, " \\\n");
    Wr.PutText(wr, "  ");
    Wr.PutText(wr, key);
    Wr.PutText(wr, " ");
    Wr.PutText(wr, Fmt.Int(arg));
  END PrintIntArg;

PROCEDURE DoPrintTextArg(wr: Wr.T; t: TEXT) =
  CONST 
    ShellHotChars = SET OF CHAR{'\'', '\\', '!', '\n', '\r'};
    PlainChars = SET OF CHAR{
      'A'..'Z', '0'..'9', 'a'..'z', 
      '@', '/', ':', '_', '-', '=', '+', '.'
    };
  VAR plain: BOOL := TRUE;
  <* FATAL Thread.Alerted, Wr.Failure *>
  BEGIN
    WITH
      n = Text.Length(t)
    DO
      FOR i := 0 TO n-1 DO 
        plain := plain AND Text.GetChar(t, i) IN PlainChars;
      END;
      IF plain THEN 
        Wr.PutText(wr, t)
      ELSE
        Wr.PutChar(wr, '\'');
        FOR i := 0 TO n-1 DO 
          WITH c = Text.GetChar(t, i) DO
            IF c IN ShellHotChars THEN
              Wr.PutChar(wr, '\\');
              Wr.PutChar(wr, c)
            ELSE
              Wr.PutChar(wr, c)
            END
          END
        END;
        Wr.PutChar(wr, '\'');
      END
    END
  END DoPrintTextArg;

BEGIN
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
