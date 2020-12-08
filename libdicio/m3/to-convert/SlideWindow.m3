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

MODULE SlideWindow EXPORTS Main;

IMPORT Rd, Wr, Thread, Fmt, Text, Char, Scan, Process;
IMPORT ParseParams AS PP;
IMPORT Basics, Reduced, Encoding, PlainEncoding;
FROM Basics IMPORT NAT, POS;
FROM Reduced IMPORT Letter;
FROM Stdio IMPORT stdin, stdout, stderr;

CONST MaxWidth = 100;

TYPE 
  Options = RECORD
    width: POS;                (* Width of sliding window *)
    commentChars: SET OF CHAR; (* Start-of-comment characters *)
    validChars: SET OF CHAR;   (* Valid characters *)
    ignoreChars: SET OF CHAR;  (* Characters to ignore *)
  END;

EXCEPTION
  MissingFinalEOL;

VAR o: Options;

PROCEDURE Main() =

  VAR encoding: PlainEncoding.T := PlainEncoding.New();

      nLettersIn: CARDINAL := 0;
      
  PROCEDURE NextLetter(): Letter RAISES {Rd.EndOfFile} =
    VAR c: CHAR;
    BEGIN
      TRY
        c := NextSignifChar(stdin);
        IF NOT c IN o.validChars THEN 
          ErrorBadChar(c, "stdin", stdin);
          <* ASSERT FALSE *>
        END;
        INC(nLettersIn);
        RETURN encoding.CharToLetter(c);
      EXCEPT
      | MissingFinalEOL => 
          ErrorMissingFinalEOL("stdin");
          <* ASSERT FALSE *>
      | Encoding.BadChars, Encoding.BadChar => 
          ErrorBadChar(c, "stdin", stdin);
          <* ASSERT FALSE *>
      END;
    END NextLetter;
    
  VAR pos: NAT;  (* buf[pos..pos+width-1] is next word. Indices are mod $width$ *)

  VAR nLettersOut: CARDINAL := 0;   (* Num letters written to stdout *)
      nWordsOut: CARDINAL := 0;     (* Num words written to stdout *)

  <* FATAL Rd.Failure, Wr.Failure, Thread.Alerted *>
  <* FATAL Encoding.BadChar, Encoding.BadLetter *>
  BEGIN
    o := GetOptions();
    (* Get pattern automaton: *)

    WITH 
      buf = NEW(REF ARRAY OF Letter, o.width)^
    DO
      TRY
        (* Initialize buffer: *)
        FOR i := 0 TO o.width-1 DO 
          buf[i] := NextLetter();
        END;
        pos := 0;

        LOOP
          (* Output next word: *)
          VAR j: NAT := pos;
          BEGIN
            FOR i := 0 TO o.width-1 DO
              Wr.PutChar(stdout, encoding.LetterToChar(buf[j]));
              INC(j); IF j >= o.width THEN j := 0 END;
            END;
          END;
          Wr.PutChar(stdout, '\n');
          INC(nWordsOut);
          INC(nLettersOut, o.width);
          
          (* Read next character: *)
          buf[pos] := NextLetter();
          INC(pos);  IF pos >= o.width THEN pos := 0 END;
        END;
      EXCEPT
      Rd.EndOfFile => (* Ok *)
      END;
    END;
    
    (* Flush and close files: *)
    Rd.Close(stdin);
    Wr.Close(stdout);
    
    (* Print statistics *)
    Wr.PutText(stderr, Fmt.Pad(Fmt.Int(nLettersIn), 10)  & " letters read\n");
    Wr.PutText(stderr, Fmt.Pad(Fmt.Int(nLettersOut), 10) & " letters written\n");
    Wr.PutText(stderr, Fmt.Pad(Fmt.Int(nWordsOut), 10)   & " words written\n");
      
    Wr.Flush(stderr);
  END Main;
  
PROCEDURE NextSignifChar (rd: Rd.T): CHAR RAISES {Rd.EndOfFile, MissingFinalEOL} =
  (* Returns next character from "rd", skipping comments and fillers. *)
  VAR c: CHAR;
      afterEOL: BOOLEAN := FALSE;
  <* FATAL Rd.Failure, Thread.Alerted *>
  BEGIN
    LOOP
      TRY
        c := Rd.GetChar(rd);
      EXCEPT 
      | Rd.EndOfFile =>
          IF afterEOL OR (Rd.Index(rd) = 0) THEN
            RAISE Rd.EndOfFile
          ELSE
            RAISE MissingFinalEOL
          END;
      END;
      
      IF c = '\n' THEN afterEOL := TRUE END;

      IF c IN o.commentChars THEN
        (* Skip comment text until newline *)
        REPEAT 
          TRY
            c := Rd.GetChar(rd);
          EXCEPT
            Rd.EndOfFile => RAISE MissingFinalEOL
          END;
        UNTIL c = '\n';
      ELSIF c IN o.ignoreChars THEN
        (* Ignore *)
      ELSE
        RETURN c
      END
    END
  END NextSignifChar;

PROCEDURE GetOptions(): Options =

  PROCEDURE GetNextCharSet(): SET OF CHAR RAISES {Scan.BadFormat} =
    VAR chars: SET OF CHAR := SET OF CHAR{};
    BEGIN
      WITH arg = PP.GetNext() DO
        FOR i := 0 TO Text.Length(arg) - 1 DO
          WITH c = Text.GetChar(arg, i) DO
            chars := chars + SET OF CHAR{c};
          END
        END
      END;
      RETURN chars
    END GetNextCharSet;
    
  VAR o: Options;

  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    TRY
      PP.BeginParsing(stderr);
      
      PP.GetKeyword("-width");
      o.width := PP.GetNextInt(1,MaxWidth);

      o.ignoreChars := Char.Spaces;
      o.commentChars := SET OF CHAR {'#'};
      o.ignoreChars := Char.Spaces;

      o.validChars := Char.Graphics + SET OF CHAR{'\200'..'\377'};
      o.validChars := o.validChars - o.commentChars - o.ignoreChars;

      IF PP.KeywordPresent("-commentChars") THEN 
        o.commentChars := GetNextCharSet();
      END;
      
      IF PP.KeywordPresent("-ignoreChars") THEN 
        o.ignoreChars := GetNextCharSet();
      END;

      IF PP.KeywordPresent("-validChars") THEN 
        o.validChars := GetNextCharSet();
      END;

      IF o.validChars = SET OF CHAR{} THEN
        Wr.PutText(stderr, "GetOptions: -validChars is empty!\n");
        RAISE Scan.BadFormat
      END;

      IF o.ignoreChars * o.validChars # SET OF CHAR {}
      OR o.ignoreChars * o.commentChars # SET OF CHAR {}
      OR o.validChars * o.commentChars # SET OF CHAR{} THEN
        Wr.PutText(stderr, "GetOptions: conflicting character set specs\n");
        RAISE Scan.BadFormat
      END;

      IF '\n' IN o.validChars OR '\n' IN o.commentChars THEN
        Wr.PutText(stderr, 
          "GetOptions: -validChars or -commentChars should not include <newline>.\n"
        );
        RAISE Scan.BadFormat
      END;
      o.ignoreChars := o.ignoreChars + SET OF CHAR{'\n'};

      PP.UnparsedTail();
      PP.EndParsing()
    EXCEPT
    | ParseParams.Error => 
        Wr.PutText(stderr, "usage: \n");
        Wr.PutText(stderr, "  slidewindow \\\n");
        Wr.PutText(stderr, "    -width n \\\n");
        Wr.PutText(stderr, "    [-commentChars chars] \\\n");
        Wr.PutText(stderr, "    [-validChars chars]\\\n");
        Wr.PutText(stderr, "    [-ignoreChars chars]\\\n");
        Wr.PutText(stderr, "      < infile > outfile\\\n");
        Wr.Flush(stderr);
        Process.Exit(1)
    END;
    RETURN o
  END GetOptions;
  
PROCEDURE ErrorMissingFinalEOL(file: TEXT) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(stderr, 
      "missing final newline in file \"" & file & "\"\n"
    );
    Wr.Flush(stderr);
    Process.Exit(1)
  END ErrorMissingFinalEOL;

PROCEDURE ErrorBadChar(c: CHAR; file: TEXT; rd: Rd.T) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(stderr, "Invalid character \'");
    IF c IN Char.Graphics THEN
      Wr.PutChar(stderr, c)
    ELSE
      Wr.PutChar(stderr, '\\');
      Wr.PutText(stderr, Fmt.Pad(Fmt.Int(ORD(c), base := 8), 3, '0'));
    END;
    Wr.PutText(stderr, "\' in file \"");
    Wr.PutText(stderr, file);
    Wr.PutText(stderr, "\" at position ");
    Wr.PutText(stderr, Fmt.Int(Rd.Index(rd)));
    Wr.PutText(stderr, "\n");
    Wr.Flush(stderr);
    Process.Exit(1)
  END ErrorBadChar;

BEGIN
  Main();
END SlideWindow.
