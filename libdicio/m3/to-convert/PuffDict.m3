(****************************************************************************)
(* (C) Copyright 1992 Universidade Estadual de Campinas (UNICAMP)           *)
(*                    Campinas, SP, Brazil                                  *)
(*                                                                          *)
(* Authors:                                                                 *)
(*                                                                          *)
(*   Tomasz Kowaltowski  - CS Dept, UNICAMP <tomasz@dcc.unicamp.ansp.br>    *)
(*   Claudio L. Lucchesi - CS Dept, UNICAMP <lucchesi@dcc.unicamp.ansp.br>  *)
(*   Jorge Stolfi        - DEC Systems Research Center <stolfi@src.dec.com> *)
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

(* Last modified on Thu Mar 19 12:38:27 PST 1992 by stolfi                  *)

MODULE PuffDict EXPORTS Main;

(*
  This program reads a word list that was compressed with
  ShrinkDict, and expands it back to its original form.
  *)

IMPORT Rd, Wr, Fmt, Thread;
FROM Stdio IMPORT stdin, stdout, stderr;

CONST MaxPrefixLength = 36; (* 10 digits + 26 letters *)

PROCEDURE Main() =
  VAR s: ARRAY [0..MaxPrefixLength-1] OF CHAR;
      s_len: [0..MaxPrefixLength] := 0;
      c: CHAR;
      n: CARDINAL;
      n_words: CARDINAL := 0;
      
  <* FATAL Wr.Failure, Rd.Failure, Thread.Alerted *>
  BEGIN
    TRY
      LOOP
        (* Here s[0..s_len-1] is the previous expanded word, *)
        (*   truncated to MaxPrefixLength bytes.             *)
        
        (* Read first character of compressed word: *)
        c := Rd.GetChar(stdin);
        INC(n_words);

        (* Decode c into prefix length n: *)
        IF c >= '0' AND c <= '9' THEN
          n := ORD(c) - ORD('0')
        ELSIF c >= 'A' AND c <= 'Z' THEN
          n := ORD(c) - ORD('A') + 10
        ELSE
          Wr.PutText(stderr, 
            "Invalid prefix length code at byte " & 
            Fmt.Int(Rd.Index(stdin)) &
            " = \'\\" & Fmt.Pad(Fmt.Int(ORD(c), base := 8), 3, '0') &
            "\'\n"
          );
          Wr.Flush(stderr);
          <* ASSERT FALSE *>
        END;

        IF n > s_len THEN
          Wr.PutText(stderr, 
            "Inconsistent prefix length at byte " & 
            Fmt.Int(Rd.Index(stdin)) &
            " = \'\\" & Fmt.Pad(Fmt.Int(ORD(c), base := 8), 3, '0') &
            "\'\n"
          );
          Wr.Flush(stderr);
          <* ASSERT FALSE *>
        END;
        
        (* Copy first n characters of previous word to output: *)
        FOR i := 0 TO n-1 DO Wr.PutChar(stdout, s[i]) END;
        
        (* Now copy rest of line and save in s[]: *)
        LOOP
          TRY c := Rd.GetChar(stdin) EXCEPT Rd.EndOfFile => c := '\n' END;
          Wr.PutChar(stdout, c);
          IF c = '\n' THEN EXIT END;
          IF n <= LAST(s) THEN s[n] := c; INC(n) END;
        END;
        (* Update string length: *)
        s_len := n;
      END
    EXCEPT
    | Rd.EndOfFile => (* Ok *)
    END;
    Wr.Close(stdout);
    Wr.PutText(stderr, Fmt.Int(n_words) & " words read\n");
  END Main;

BEGIN
  Main()
END PuffDict.
