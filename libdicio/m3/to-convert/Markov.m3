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

(* Last modified on Thu Aug 20 17:01:51 PDT 1992 by stolfi                  *)

MODULE Markov EXPORTS Main;

IMPORT
  Rd, Wr, Text, Fmt, Char, Random, FileStream, ParseParams,
  Scan, RTMisc, Thread;
FROM Stdio IMPORT stdin, stderr, stdout;

TYPE
  NAT = CARDINAL;
  BOOL = BOOLEAN;

  Options = RECORD
      order: [0..2];         (* Number of elements in tuple *)
      positional: BOOL;      (* TRUE to compute separate tuple tables by column *)
      nWords: NAT;           (* Number of words to generate *)
      tableFileName: TEXT;   (* File where to write the table, or "" if stderr *)
    END;

CONST
  MaxChars = 63;       (* Max distinct chars in input *)
  MaxWordLength = 40;  (* Maximum input word length *)

TYPE
  Code = [0..MaxChars];       (* Code 0 is EOL (filler at begin/end of word) *)

  Col = [0..MaxWordLength];   (* Column 0 is first letter, Length(word) is the EOL *)

  CharMap = REF RECORD
      nChars: NAT;             (* Number of distinct chars seen in input *)
      dir: ARRAY CHAR OF Code; (* Internal Code of each input CHAR *)
      inv: ARRAY Code OF CHAR; (* External representation of each Code *)
    END;

TYPE
  Range = RECORD colLast, aLast, bLast, cLast: NAT END;

  Page = REF ARRAY Code OF NAT; (* Frequency of each character code *)

  PageIndex = NAT;              (* Index of a tuple minus last component *)

  Table = REF RECORD
      page: REF ARRAY OF Page;  (* Maps tuple page to CodeFreqTable *)
      map: CharMap;             (* Mapping of CHAR to Code and vice versa *)
      range: Range;             (* Range of tuples actually present in table *)
      nPages: NAT;              (* Number of pages (contexts) seen *)
      nTuples: NAT;             (* Number of distinct tuples seen *)
      order: NAT (*[0..2]*);    (* Number of elements in tuple *)
      positional: BOOL;         (* TRUE to add column index in front of every tuple *)
    END;

EXCEPTION
  TooManyChars; (* Too many distinct characters in input *)
  WordTooLong;  (* Word has more than MaxChars characters (positional model only) *)
  Done;         (* End of file, etc. *)

VAR o: Options;

PROCEDURE Main() =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    o := GetOptions();
    WITH
      t = CountTuples(stdin, o.order, o.positional)
    DO

      WITH
        twr = OpenWrite(o.tableFileName, def := stderr)
      DO
        PrintTable(twr, t);
      END;

      IF o.nWords > 0 THEN
        GenerateSample(stdout, t, o.nWords);
      END;

    END;
    Wr.PutText(stderr, "done.\n");
    Wr.Flush(stderr);
    Wr.Flush(stdout);
  END Main;

PROCEDURE GetOptions(): Options =
  VAR o: Options;
  <* FATAL Scan.BadFormat, Wr.Failure, Thread.Alerted *>
  BEGIN
    TRY
      ParseParams.BeginParsing(stderr);

        IF ParseParams.KeywordPresent("-table") THEN
          o.tableFileName := ParseParams.GetNext()
        ELSE
          o.tableFileName := ""
        END;

        ParseParams.GetKeyword("-order");
        o.order := ParseParams.GetNextInt(0, 2);

        o.positional := ParseParams.KeywordPresent("-positional");

        IF ParseParams.KeywordPresent("-nWords") THEN
          o.nWords := ParseParams.GetNextInt(0, 1000000)
        ELSE
          o.nWords := 0
        END;

      ParseParams.EndParsing();
    EXCEPT
    | Scan.BadFormat =>
        Wr.PutText(stderr, "Usage: \\\n");
        Wr.PutText(stderr, "  Markov \\\n");
        Wr.PutText(stderr, "    -order { 0 | 1 | 2 } \\\n");
        Wr.PutText(stderr, "    [ -positional ] \\\n");
        Wr.PutText(stderr, "    [ -nWords <nnn> ] \\\n");
        Wr.PutText(stderr, "    [ -table <table-file-name> ] \\\n");
        Wr.PutText(stderr, "  < original.dic  \\\n");
        Wr.PutText(stderr, "  > synthetic.dic \\\n");
        RTMisc.Exit(1);
    END;
    RETURN o
  END GetOptions;

PROCEDURE OpenWrite(name: TEXT; def: Wr.T := NIL): Wr.T =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    IF Text.Empty(name) THEN
      RETURN def
    ELSE
      RETURN FileStream.OpenWrite(name)
    END;
  END OpenWrite;

PROCEDURE GetPage(t: Table; col: NAT; a, b: Code; create: BOOL := FALSE): Page =
(*
  Returns the page of the table that corresponds to the given context
  (column and two previous character codes).  Creates and clears the
  page if missing and /create=TRUE/. *)

  VAR i: PageIndex;
  CONST NCodes = LAST(Code) - FIRST(Code) + 1;
  BEGIN
    WITH r = t.range DO
      IF t.positional THEN <* ASSERT col <= r.colLast *> i := col ELSE i := 0 END;
      IF t.order >= 2 THEN <* ASSERT a <= r.aLast *> i := i * NCodes + a END;
      IF t.order >= 1 THEN <* ASSERT b <= r.bLast *> i := i * NCodes + b END;
    END;
    WITH page = t.page[i] DO
      IF page = NIL AND create THEN
        INC(t.nPages);
        page := NEW(Page);
        WITH f = page^ DO
          FOR c := FIRST(f) TO LAST(f) DO f[c] := 0 END
        END;
      END;
      RETURN page
    END
  END GetPage;

PROCEDURE NumPages(order: [0..2]; positional: BOOL): NAT =
  VAR n: NAT;
  CONST NCodes = LAST(Code) - FIRST(Code) + 1;
  CONST NCols = LAST(Col) - FIRST(Col) + 1;
  BEGIN
    IF positional THEN n := NCols ELSE n := 1 END;
    CASE order OF
    | 0 => RETURN n
    | 1 => RETURN n*NCodes
    | 2 => RETURN n*NCodes*NCodes
    END
  END NumPages;

PROCEDURE NewCharMap(): CharMap =
  BEGIN
    WITH
      map = NEW(CharMap),
      m = map^
    DO
      (* Initialize character map: *)
      m.nChars := 0;
      FOR ch := FIRST(CHAR) TO LAST(CHAR) DO m.dir[ch] := 0 END;
      FOR c := FIRST(Code) TO LAST(Code) DO m.inv[c] := Char.NUL END;
      RETURN map;
    END
  END NewCharMap;

PROCEDURE CharToCode(map: CharMap; ch: CHAR): Code =
  <* FATAL Wr.Failure, Thread.Alerted, TooManyChars *>
  BEGIN
    WITH m = map^, dir = m.dir, inv = m.inv DO
      (*
      Wr.PutText(stderr, "(ch = " & Fmt.Int(ORD(ch)) & ")\n");
      Wr.Flush(stderr);
      Wr.PutText(stderr, "(m.nChars = " & Fmt.Int(m.nChars) & ")\n");
      Wr.Flush(stderr);
      *)
      IF dir[ch] = 0 THEN
        IF m.nChars = MaxChars THEN RAISE TooManyChars END;
        INC(m.nChars);
        dir[ch] := m.nChars;
        inv[m.nChars] := ch;
      END;
      RETURN dir[ch]
    END;
  END CharToCode;

PROCEDURE CountTuples(rd: Rd.T; order: [0..2]; positional: BOOL): Table =
(*
  Returns a table with absolute occurence counts for character tuples of
  the given /order/ in words read from /rd/ (one word for line).
  *)

  BEGIN
    WITH
      t = NEW(Table,
        map := NewCharMap(),
        range := Range{colLast := 0, aLast := 0, bLast := 0, cLast := 0},
        nPages := 0,
        nTuples := 0,
        order := order,
        positional := positional,
        page := NEW(REF ARRAY OF Page, NumPages(order, positional))
      )
    DO

      VAR a, b, c: Code;
          col: NAT := 0;

      PROCEDURE Tally(code: Code) =
        BEGIN
          a := b; b := c; c := code;

          (* Update range of tuples: *)
          WITH r = t.range DO
            IF positional THEN r.colLast := MAX(r.colLast, col) END;
            IF order >= 2 THEN r.aLast := MAX(r.aLast , a) END;
            IF order >= 1 THEN r.bLast := MAX(r.bLast, b) END;
            r.cLast := MAX(r.cLast, c)
          END;

          (* Tally it: *)
          WITH
            page = GetPage(t, col, a, b, create := TRUE),
            fi = page^[c]
          DO
            IF fi = 0 THEN INC(t.nTuples) END;
            fi := fi + 1;
          END;

        END Tally;

      VAR ch: CHAR;
          nWords: NAT := 0;

      <* FATAL Wr.Failure, Rd.Failure, Rd.EndOfFile, Thread.Alerted *>
      <* FATAL WordTooLong *>
      BEGIN
        (* Clear frequency counts: *)
        FOR i := 0 TO LAST(t.page^) DO t.page[i] := NIL END;

        (* Count tuples: *)
        TRY
          LOOP (* On words *)
            TRY ch := Rd.GetChar(rd) EXCEPT Rd.EndOfFile => RAISE Done END;
            IF ch # '\n' THEN
              INC(nWords);
              a := 0; b := 0; c := 0;
              col := 0;
              REPEAT
                Tally(CharToCode(t.map, ch));
                ch := Rd.GetChar(rd);
                INC(col);
                IF col > MaxWordLength THEN RAISE WordTooLong END;
              UNTIL ch = '\n';
              Tally(0);
              IF nWords MOD 1000 = 0 THEN
                Wr.PutChar(stderr, '.'); Wr.Flush(stderr);
              END;
            END;
          END;

        EXCEPT
        | Done => (* OK *)
        END;

        Wr.PutChar(stderr, '\n');
        Wr.PutText(stderr, Fmt.Pad(Fmt.Int(nWords), 8) & " words read\n");
        Wr.PutText(stderr, Fmt.Pad(Fmt.Int(t.nPages), 8) & " pages in table\n");
        Wr.PutText(stderr, Fmt.Pad(Fmt.Int(t.nTuples), 8) & " distinct tuples\n");
        Wr.Flush(stderr);
      END;

      RETURN t
    END;
  END CountTuples;

PROCEDURE PrintTuple(
    wr: Wr.T;
    map: CharMap;
    order: [0..2];
    positional: BOOL;
    col: Col;
    a, b, c: Code;
  ) =

  PROCEDURE PrintChar(code: Code) =
    <* FATAL Wr.Failure, Thread.Alerted *>
    BEGIN
      IF code = 0 THEN
        Wr.PutText(wr, "EOL")
      ELSE
        WITH ch = map.inv[code] DO
          IF ch IN Char.Graphics THEN
            Wr.PutChar(wr, '`');
            Wr.PutChar(wr, ch);
            Wr.PutChar(wr, '\'')
          ELSE
            Wr.PutText(wr, Fmt.Pad(Fmt.Int(ORD(ch), base := 8), 3, '0'))
          END
        END
      END;
      Wr.PutChar(wr, ' ');
    END PrintChar;

  PROCEDURE PrintColumn(col: Col) =
    <* FATAL Wr.Failure, Thread.Alerted *>
    BEGIN
      Wr.PutText(wr, Fmt.Pad(Fmt.Int(col), 3));
      Wr.PutChar(wr, ' ');
    END PrintColumn;

  BEGIN
    WITH m = map^ DO
      IF positional THEN PrintColumn(col) END;
      IF order >= 2 THEN PrintChar(a) END;
      IF order >= 1 THEN PrintChar(b) END;
      IF order >= 0 THEN PrintChar(c) END;
    END
  END PrintTuple;

PROCEDURE PrintTable(wr: Wr.T; t: Table) =
  VAR r: Range := t.range;
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    BEGIN
      FOR col := 0 TO r.colLast DO
        FOR a := 0 TO r.aLast DO
          FOR b := 0 TO r.bLast DO
            WITH
              page = GetPage(t, col, a, b)
            DO
              IF page # NIL THEN
                WITH f = page^ DO
                  VAR sum: NAT := 0;
                  BEGIN
                    (* Compute total occurrences for tuples ([col,]a,b,ANY) *)
                    FOR c := 0 TO r.cLast DO
                      WITH fi = f[c] DO sum := sum + fi END;
                    END;
                    <* ASSERT sum # 0 *>

                    (* Print absolute counts and frequs for all observed tuples *)
                    FOR c := 0 TO r.cLast DO
                      WITH
                        fi = f[c],
                        ffi = FLOAT(fi)/FLOAT(sum)
                      DO
                        IF fi # 0 THEN
                          PrintTuple(wr, t.map, t.order, t.positional, col, a, b, c);
                          Wr.PutText(wr, Fmt.Pad(Fmt.Int(fi), 8));
                          Wr.PutText(wr, Fmt.Pad(Fmt.Real(ffi, 4, Fmt.Style.Flo), 8));
                          Wr.PutChar(wr, '\n');
                        END
                      END;
                    END;
                  END;
                END;
              END;
            END;
          END;
        END;
      END;
    END;
    Wr.Flush(wr);
  END PrintTable;

CONST ProbScale = 1000000;

PROCEDURE GenerateSample(wr: Wr.T; t: Table; nWords: NAT) =
  VAR r: Range := t.range;
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    FOR col := 0 TO r.colLast DO
      FOR a := 0 TO r.aLast DO
        FOR b := 0 TO r.bLast DO
          WITH page = GetPage(t, col, a, b) DO
            IF page # NIL THEN IntegrateAndNormalize(page^) END;
          END;
        END;
      END;
    END;

    (* Generate words with same frequencies: *)
    FOR i := 1 TO nWords DO
      GenerateWord(wr, t);
      IF i MOD 1000 = 0 THEN
        Wr.PutChar(stderr, '.'); Wr.Flush(stderr);
      END;
    END;

    Wr.Flush(wr);
  END GenerateSample;

PROCEDURE IntegrateAndNormalize(VAR f: ARRAY OF NAT) =
  VAR sum: NAT := 0;
  BEGIN
    IF NUMBER(f) = 0 THEN RETURN END;

    FOR c := 0 TO LAST(f) DO
      WITH fc = f[c] DO
        sum := sum + fc;
        fc := sum
      END;
    END;

    <* ASSERT sum # 0 *>
    FOR c := 0 TO LAST(f)-1 DO
      WITH fc = f[c] DO
        fc := ROUND(FLOAT(fc) * FLOAT(ProbScale)/FLOAT(sum))
      END;
    END;
    f[LAST(f)] := ProbScale;
  END IntegrateAndNormalize;

PROCEDURE GenerateWord(wr: Wr.T; t: Table) =
  VAR a, b, c: Code := 0;
      col: NAT := 0;
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    REPEAT
      (* Choose the /c/ code, based on /a/ and /b/: *)
      a := b; b := c;
      c := GenerateCode(t, col, a, b);
      (* Output /c/ as a CHAR: *)
      IF c = 0 THEN
        Wr.PutChar(wr, '\n')
      ELSE
        Wr.PutChar(wr, t.map.inv[c])
      END;
      INC(col)
    UNTIL c = 0;
    Wr.Flush(wr);
  END GenerateWord;

PROCEDURE GenerateCode(t: Table; col: NAT; a, b: Code): Code =
  VAR coin: NAT := Random.Subrange(NIL, 0, ProbScale-1);
      lo, hi, c: Code;
  BEGIN
    WITH f = GetPage(t, col, a, b)^ DO
      (* Binary search: *)
      hi := t.range.cLast;
      IF coin >= f[hi] THEN RETURN hi END;
      lo := 0;
      WHILE lo < hi DO
        c := (lo + hi + 1) DIV 2;
        <* ASSERT c > lo *>
        <* ASSERT c <= hi *>
        IF coin <= f[c-1] THEN
          hi := c - 1
        ELSIF coin > f[c] THEN
          lo := c + 1
        ELSE
          RETURN c
        END
      END;
      RETURN lo
    END
  END GenerateCode;

BEGIN
  Main()
END Markov.



