MODULE StringPrinter;

(* See the copyright and disclaimer note at the end of this file. *)

IMPORT Wr, Text, Thread;
IMPORT Encoding;
FROM Basics IMPORT String, Abort;
FROM Basics IMPORT NAT, BOOL;

REVEAL
  T = Public BRANDED OBJECT 

        (* parameters of Init: *)
        wr: Wr.T;             (* The underlying writer *)
        encoding: Encoding.T; (* The external encoding *)
        empty: TEXT;          (* Symbol denoting the empty string *)
        emptyLength: NAT;     (* Length of same *)
        sep: TEXT;            (* String separator *)
        sepLength: NAT;       (* Length of same *)
        maxChars: NAT;        (* Maximum total characters to print *)
        etc: TEXT;            (* Sumbol denoting "and more" *)
        etcLength: NAT;       (* Length of same *)
        rightMargin: NAT;     (* Nominal right margin *)
        leftMargin: NAT;      (* Nominal left margin *)
        initialColumn: NAT;   (* Starting column *)

        (* running counts: *)
        currentColumn: NAT;   (* Current column *)
        totChars: NAT;        (* Total characters printed since Init/Reset *)

      OVERRIDES
        PutString := PutString;
        Reset := Reset;
      END;

PROCEDURE Init(
    self: T;
    wr: Wr.T;
    encoding: Encoding.T;
    empty: TEXT := "()";
    sep: TEXT := ", ";
    maxChars: NAT := LAST(NAT);
    etc: TEXT := "...";
    rightMargin: NAT := 60;
    leftMargin: NAT := 0;
    initialColumn: NAT := 0;
  ) =
  BEGIN
    self.wr := wr;
    self.encoding := encoding;
    self.empty := empty;
    self.emptyLength := Text.Length(empty);
    self.sep := sep;
    self.sepLength := Text.Length(sep);
    self.maxChars := maxChars;
    self.etc := etc;
    self.etcLength := Text.Length(etc);
    self.rightMargin := rightMargin;
    self.leftMargin := leftMargin;
    self.initialColumn := initialColumn;
    self.Reset();
  END Init;

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
  ): T =
  BEGIN
    WITH s = NEW(T) DO
      Init(s,
        wr := wr,
        encoding := encoding,
        empty := empty,
        sep := sep,
        maxChars := maxChars,
        etc := etc,
        leftMargin := leftMargin,
        rightMargin := rightMargin,
        initialColumn := initialColumn
      );
      RETURN s
    END;
  END New;

PROCEDURE Reset(self: T) =
  BEGIN
    self.currentColumn := self.initialColumn;
    self.totChars := 0;
  END Reset;

PROCEDURE PutString(self: T; READONLY w: String; rev: BOOL := FALSE) RAISES {Abort} =
  VAR wLength: NAT;
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    IF self.totChars # 0 THEN
      Wr.PutText(self.wr, self.sep);
      INC(self.totChars, self.sepLength);
      INC(self.currentColumn, self.sepLength);
    END;
    IF NUMBER(w) = 0 THEN
      wLength := self.emptyLength
    ELSE
      wLength := NUMBER(w)
    END;
    IF self.totChars + wLength > self.maxChars THEN
      Wr.PutText(self.wr, self.etc);
      INC(self.totChars, self.etcLength);
      INC(self.currentColumn, self.etcLength);
      RAISE Abort
    END;
    IF self.currentColumn + wLength > self.rightMargin THEN
      Wr.PutChar(self.wr, '\n');
      FOR i := 0 TO self.leftMargin - 1 DO 
        Wr.PutChar(self.wr, ' ')
      END;
      self.currentColumn := self.leftMargin
    END;
    IF NUMBER(w) = 0 THEN
      Wr.PutText(self.wr, self.empty)
    ELSIF rev THEN
      FOR i := LAST(w) TO 0 BY -1 DO 
        self.encoding.PrintLetter(self.wr, w[i])
      END;
    ELSE
      FOR i := 0 TO LAST(w) DO
        self.encoding.PrintLetter(self.wr, w[i])
      END;
    END;
    INC(self.totChars, wLength);
    INC(self.currentColumn, wLength);
  END PutString;

BEGIN
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
