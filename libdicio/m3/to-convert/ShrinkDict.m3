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

(* Last modified on Thu Mar 19 12:38:42 PST 1992 by stolfi                  *)

MODULE ShrinkDict EXPORTS Main;

(*
  This program reads a list of words from stdin,
  preferably in alphabetic order, and writes out a version
  where repeated prefixes are coded by single characters. 
  
  Specifically, each line of the output starts with a
  count n in [0..36], encoded as one character 0-9A-Z.
  This character is a surrogate for the first n characters of the 
  preceding word.  (The first line always starts with a '0'). *)

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
    s[0] := '\n';
    TRY
      LOOP
        (* Here s[0..s_len-1] is the previous word, truncated to MaxPrefixLength bytes. *)
        (* Skip prefix of word that is shared with previous word: *)
        n := 0;
        c := Rd.GetChar(stdin);
        INC(n_words);
        WHILE (c # '\n') AND (n < s_len) AND (c = s[n]) DO
          INC(n); 
          TRY c := Rd.GetChar(stdin) EXCEPT Rd.EndOfFile => c := '\n' END
        END;
        (* Output length of prefix n, encoded as one character: *)
        IF n < 10 THEN
          Wr.PutChar(stdout, VAL(ORD('0') + n, CHAR))
        ELSE
          Wr.PutChar(stdout, VAL(ORD('A') + (n - 10), CHAR))
        END;
        (* Copy remaining characters, and save them in s[]: *)
        WHILE c # '\n' DO
          Wr.PutChar(stdout, c);
          IF n <= LAST(s) THEN s[n] := c; INC(n) END;
          TRY c := Rd.GetChar(stdin) EXCEPT Rd.EndOfFile => c := '\n' END
        END;
        Wr.PutChar(stdout, '\n');
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
END ShrinkDict.
