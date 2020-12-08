MODULE PlainEncoding;

(* See the copyright and disclaimer note at the end of this file. *)

IMPORT Rd, Wr, TextF, TextWr, Fmt, Thread; 
IMPORT Basics;

FROM Basics IMPORT Symbol, String, Done;
FROM Basics IMPORT NAT, CHARS;
FROM Encoding IMPORT BadChar, BadChars, BadLetter, BadString;

REVEAL 
  T = Public BRANDED 
    OBJECT 
    OVERRIDES
      StringToChars := StringToChars;
      StringToText := StringToText;
      PrintString := PrintString;
      PrintLetter := PrintLetter;

      CharsToString := CharsToString;
      TextToString := TextToString;
      ReadString := ReadString;

      CharToLetter := CharToLetter;
      SymbolToChar := SymbolToChar;
    END;

PROCEDURE StringToChars(
    e: T;
    READONLY w: String; 
    VAR (*IO*) c: REF CHARS; 
    VAR (*OUT*) len: NAT;
  ) RAISES {BadString} =
  BEGIN
    len := NUMBER(w);
    IF c = NIL OR NUMBER(c^) < len THEN Basics.ExpandChars(c, len) END;
    WITH cc = c^ DO
      FOR i := 0 TO LAST(w) DO 
        TRY
          cc[i] := SymbolToChar(e, w[i])
        EXCEPT
        | BadLetter => RaiseBadString(w[i], i)
        END 
      END
    END
  END StringToChars;

PROCEDURE StringToText(
    e: T;
    READONLY w: String;
  ): TEXT RAISES {BadString} =
  BEGIN
    WITH wr = TextWr.New() DO
      PrintString(e, wr, w);
      RETURN TextWr.ToText(wr)
    END
  END StringToText;

PROCEDURE PrintString(
    <*UNUSED*> e: T;
    wr: Wr.T;
    READONLY w: String;
  ) RAISES {BadString} =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    FOR i := 0 TO LAST(w) DO 
      WITH wi = w[i] DO
        IF NOT (wi IN ValidLetters) THEN RaiseBadString(wi, i) END;
        Wr.PutChar(wr, VAL(wi, CHAR))
      END
    END
  END PrintString;

PROCEDURE RaiseBadString(let: Symbol; i: NAT) RAISES {BadString} =
  BEGIN
    RAISE BadString(
      "bad symbol at position " & Fmt.Int(i) &
      " = " & Fmt.Int(let) & " = 8_" & Fmt.Int(let, base := 8)
    )
  END RaiseBadString;

CONST OctalDigit = ARRAY [0..7] OF CHAR{'0', '1', '2', '3', '4', '5', '6', '7'};

PROCEDURE PrintLetter(
    e: T;
    wr: Wr.T; 
    symbol: Symbol;
  ) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    TRY
      Wr.PutChar(wr, SymbolToChar(e, symbol))
    EXCEPT
    | BadLetter(bad) =>
        Wr.PutChar(wr, '\\');
        WITH n = ORD(bad) DO
          Wr.PutChar(wr, OctalDigit[n DIV 64]);          
          Wr.PutChar(wr, OctalDigit[(n DIV 8) MOD 8]);          
          Wr.PutChar(wr, OctalDigit[n MOD 8]);          
        END;
    END;
  END PrintLetter;

PROCEDURE CharsToString(
    e: T;
    READONLY c: CHARS;
    VAR (*IO*) w: REF String;
    VAR (*OUT*) len: NAT;
  ) RAISES {BadChars} =
  BEGIN
    len := NUMBER(c);
    IF w = NIL OR NUMBER(w^) < len THEN Basics.ExpandString(w, len) END;
    WITH ww = w^ DO
      FOR i := 0 TO LAST(c) DO 
        TRY
          ww[i] := CharToLetter(e, c[i])
        EXCEPT
        | BadChar => RaiseBadChars(c[i], i)
        END 
      END
    END
  END CharsToString;

PROCEDURE TextToString(
    e: T;
    t: TEXT;
    VAR (*IO*) w: REF String;
    VAR (*OUT*) len: NAT;
  ) RAISES {BadChars} =
  BEGIN
    WITH tt = t^ DO
      <* ASSERT tt[LAST(tt)] = '\000' *>
      CharsToString(e, SUBARRAY(tt, 0, NUMBER(tt)-1), w, len)
    END;
  END TextToString;
  
CONST SpaceLetter: Symbol = ORD(' ');

PROCEDURE ReadString(
    e: T;
    rd: Rd.T;
    VAR (*IO*) w: REF String;
    VAR (*OUT*) len: NAT;
  ) RAISES {BadChars, Done, Rd.Failure, Thread.Alerted} =
  VAR c: CHAR;
  BEGIN
    len := 0;
    TRY 
      REPEAT c := Rd.GetChar(rd) UNTIL c # ' '
    EXCEPT 
      Rd.EndOfFile => RAISE Done 
    END;
    WHILE c # '\n' DO
      Basics.ExpandString(w, len+1);
      TRY
        w[len] := CharToLetter(e, c)
      EXCEPT 
      | BadChar => RaiseBadChars(c, len)
      END;
      INC(len);
      TRY 
        c := Rd.GetChar(rd)
      EXCEPT 
        Rd.EndOfFile => c := '\n'
      END;
    END;
    WHILE len > 0 AND w[len-1] = SpaceLetter DO DEC(len) END;
  END ReadString;

PROCEDURE RaiseBadChars(c: CHAR; i: NAT) RAISES {BadChars} =
  BEGIN
    RAISE BadChars(
      "bad character at position " & Fmt.Int(i) &
      " = \\" & Fmt.Int(ORD(c), base := 8)
    )
  END RaiseBadChars;

PROCEDURE CharToLetter(
    <*UNUSED*> e: T; 
    ch: CHAR;
  ): Symbol RAISES {BadChar} =
  BEGIN
    IF NOT(ch IN ValidChars) THEN RAISE BadChar(ch) END;
    RETURN ORD(ch)
  END CharToLetter;

PROCEDURE SymbolToChar(
    <*UNUSED*> e: T; 
    let: Symbol;
  ): CHAR RAISES {BadLetter} =
  BEGIN
    IF NOT(let IN ValidLetters) THEN RAISE BadLetter(let) END;
    RETURN VAL(let, CHAR)
  END SymbolToChar;

PROCEDURE Init(<*UNUSED*> e: T) =
  BEGIN
  END Init;

PROCEDURE New(): T =
  BEGIN
    WITH e = NEW(T) DO 
      Init(e);
      RETURN e
    END;
  END New;

BEGIN
END PlainEncoding.

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
