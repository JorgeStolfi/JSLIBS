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

MODULE PatSplit EXPORTS Main;

IMPORT Rd, Wr, Thread, FileRd, FileWr, OSError, Fmt, Text, Char, Scan, Process;
IMPORT ParseParams AS PP;
IMPORT Basics, Reduced, Encoding, PlainEncoding;
FROM Basics IMPORT NAT, BOOL, Done, Abort;
FROM Reduced IMPORT Letter, String, State, NullState, UnitState;
FROM Stdio IMPORT stdin, stdout, stderr;

TYPE 
  Options = RECORD
    commentChars: SET OF CHAR; (* Start-of-comment characters *)
    validChars: SET OF CHAR;   (* Valid characters *)
    ignoreChars: SET OF CHAR;  (* Characters to ignore *)
    breakChar: CHAR;           (* Character used in automaton to indicate a break *)
    patterns: TEXT;            (* File containing patterns, or "" *)
    patAutomaton: TEXT;        (* Dumpfile containing pattern automaton, or "" *)
    verbose: BOOL;             (* TRUE to print matched patterns *)
  END;

EXCEPTION
  MissingFinalEOL;

VAR o: Options;

PROCEDURE Main() =

  VAR encoding: PlainEncoding.T := PlainEncoding.New();

      nLettersIn: CARDINAL := 0;
      exhausted: BOOL := FALSE; (* TRUE if input file is exhausted *)

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
    
  VAR pat: Reduced.T;
      breakLetter: Letter;
  
      bufSize: CARDINAL := 256;  (* same as NUMBER(buf^) *)
      bufCount: CARDINAL := 0;   (* Num of valid letters in "buf" *)
    
      buf: REF ARRAY OF Letter := NEW(REF ARRAY OF Letter, bufSize);
        (* buf[0..bufCount-1] are the letters already read and not written *)

      mark:  REF ARRAY OF BOOL := NEW(REF ARRAY OF BOOL, bufSize);
        (* If mark[i] is TRUE, must output a newline before buf[i] *)
        
      bcount: REF ARRAY OF NAT := NEW(REF ARRAY OF NAT, bufSize);
        (* During matching, bcount[i] is the number of "breakLetter" *)
        (*   transitions taken just before reading buf[i], *)
        (*   along the current path in the pattern automaton. *)
        
  PROCEDURE MatchAndMark(i: CARDINAL; s: State; bc: NAT;): BOOL =
    (* Enumerates all matches of buf[i..] against the *)
    (*   sub-automaton rooted at "s", setting "mark" as appropriate. *)
    (* "bc" is the number of "breakLetter" transitions taken since *)
    (*   gobbling previous input letter. *)
    (* Returns TRUE if there is any match, FALSE otherwise. *)
    VAR matched: BOOL;
    BEGIN
      IF s = NullState THEN
        RETURN FALSE
      ELSIF s = UnitState THEN
        IF o.verbose THEN PrintPath(i, bc) END;
        RETURN TRUE
      ELSE
        (* Count the empty match, if any: *)
        matched := pat.Final(s);
        IF matched AND o.verbose THEN PrintPath(i, bc) END;

        (* Ensure that buf[i], mark[i], and bcount[i] exist: *)
        IF i >= bufSize THEN
          (* Expand "buf, mark": *)
          WITH 
            nbufSize = 2 * bufSize,
            nbuf = NEW(REF ARRAY OF Letter, nbufSize),
            nmark = NEW(REF ARRAY OF BOOL, nbufSize),
            nbcount = NEW(REF ARRAY OF NAT, nbufSize)
          DO
            SUBARRAY(nbuf^, 0, bufSize) := buf^;    buf := nbuf;
            SUBARRAY(nmark^, 0, bufSize) := mark^; mark := nmark;
            SUBARRAY(nbcount^, 0, bufSize) := bcount^; bcount := nbcount;
            bufSize := nbufSize;
          END;
        END;
        
        (* Ensure that buf[i] is read, or set "exhausted": *)
        IF i = bufCount AND NOT exhausted THEN
          mark[bufCount] := FALSE;
          TRY
            buf[bufCount] := NextLetter();
            INC(bufCount);
          EXCEPT
          | Rd.EndOfFile => exhausted := TRUE;
          END;
        END;
        
        (* Enumerate matches that start with break: *)
        WITH 
          t = pat.Succ(s, breakLetter),
          breakMatched = MatchAndMark (i, t, bc + 1)
        DO
          IF breakMatched THEN
            matched := TRUE;
            mark[i] := TRUE
          END
        END;
        
        (* Enumerate matches that start with next letter: *)
        IF i < bufCount THEN
          bcount[i] := bc;
          WITH 
            t = pat.Succ(s, buf[i]),
            restMatched = MatchAndMark (i+1, t, 0)
          DO
            IF restMatched THEN
              matched := TRUE;
            END
          END;
        END;
        
        (* Done: *)
        RETURN matched
      END
    END MatchAndMark;
    
  PROCEDURE PrintPath(n: CARDINAL; bc: CARDINAL) =
    (* Prints the current path in the pattern automaton. *)
    (* The path consists of the letters buf[0..n-1]; *)
    (*   each letter buf[i] is preceded by bcount[i] break marks, *)
    (*   and the last letter is followed by "bc" break marks. *)
    <* FATAL Encoding.BadLetter *>
    BEGIN
      Wr.PutText(stderr, "**  ");
      Wr.PutText(stderr, Fmt.Pad(Fmt.Int(nLettersIn - bufCount), 10));
      Wr.PutChar(stderr, ' ');
      FOR i := 0 TO n-1 DO
        WITH bci = bcount[i] DO
          FOR j := 0 TO bci - 1 DO
            Wr.PutChar(stderr, o.breakChar)
          END;
        END;
        Wr.PutChar(stderr, encoding.LetterToChar(buf[i]));
      END;
      FOR j := 0 TO bc - 1 DO
        Wr.PutChar(stderr, o.breakChar)
      END;
      Wr.PutChar(stderr, '\n');
    END PrintPath;

  VAR nLettersOut: CARDINAL := 0;   (* Num letters written to stdout *)
      nWordsOut: CARDINAL := 0;     (* Num words written to stdout *)
      curWordLength: CARDINAL := 0; (* Letters written since last newline *)
      root: State;

  <* FATAL Rd.Failure, Wr.Failure, Thread.Alerted *>
  <* FATAL Encoding.BadChar, Encoding.BadLetter *>
  BEGIN
    o := GetOptions();
    (* Get pattern automaton: *)
    pat := GetAutomaton(o.patterns, o.patAutomaton);
    root := pat.Root();
    breakLetter := encoding.CharToLetter(o.breakChar);

    WHILE bufCount > 0 OR NOT exhausted DO
      (* Match patterns against remaining stdin chars,         *)
      (*   starting at current pos, and mark all breaks found: *)
      EVAL MatchAndMark(0, root, 0);

      (* Output a newline, if this letter got marked: *)
      IF mark[0] AND curWordLength > 0 THEN
        Wr.PutChar(stdout, '\n');
        INC(nWordsOut);
        curWordLength := 0;
        Wr.Flush(stdout);
      END;
      
      (* Output this letter: *)
      Wr.PutChar(stdout, encoding.LetterToChar(buf[0]));
      INC(nLettersOut);
      INC(curWordLength);
      
      (* Discard first letter and shift "buf, mark": *)
      FOR i := 0 TO bufCount - 2 DO
        buf[i] := buf[i+1];
        mark[i] := mark[i+1]
      END;
      DEC(bufCount);
    END;
    
    (* Write last word, if any: *)
    IF curWordLength > 0 THEN
      Wr.PutChar(stdout, '\n');
      INC(nWordsOut);
      Wr.Flush(stdout);
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
  
PROCEDURE GetAutomaton(patFile, dmpFile: TEXT): Reduced.T =
  BEGIN
    IF NOT Text.Empty(patFile) THEN
      
      VAR e: Encoding.T := PlainEncoding.New();
          rd: Rd.T := OpenFile(patFile);
      
      PROCEDURE NextPat (* : Reduced.NextStringProc *) (
          VAR (*IO*) s: REF String;  (* The string buffer *)
          VAR (*OUT*) len: NAT;      (* Length of returned string *)
          VAR (*OUT*) add: BOOLEAN;  (* TRUE to add, FALSE to delete *)
        ) RAISES {Done, Abort} =
        <* FATAL Encoding.BadChars, Rd.Failure, Thread.Alerted *>
        BEGIN
          e.ReadString(rd, s, len);
          add := TRUE
        END NextPat;
        
      <* FATAL Basics.Abort *>
      BEGIN
        WITH aut = Reduced.New(size := 1000) DO
          aut.Build(next := NextPat);
          RETURN aut
        END
      END
    ELSIF NOT Text.Empty(dmpFile) THEN
      RETURN Reduced.Load(OpenFile(dmpFile))
    ELSE
      <* ASSERT FALSE *>
    END;
  END GetAutomaton;

PROCEDURE OpenFile(file: TEXT): Rd.T =
  BEGIN
    IF Text.Equal(file, "-") THEN
      RETURN stdin
    ELSE
      TRY
        RETURN FileRd.Open(file);
      EXCEPT
      | Rd.Failure => ErrorCannotOpen(file); <* ASSERT FALSE *>
      END
    END
  END OpenFile;
  
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

      IF PP.KeywordPresent("-patterns") THEN
        o.patterns := PP.GetNext();
      ELSE
        o.patterns := "";
      END;
      IF PP.KeywordPresent("-patAutomaton") THEN
        o.patAutomaton := PP.GetNext();
      ELSE
        o.patAutomaton := "";
      END;
      
      IF Text.Empty(o.patterns) = Text.Empty(o.patAutomaton) THEN
        Wr.PutText(stderr, "GetOptions: exactly one of \"-patterns\"");
        Wr.PutText(stderr, " and \"-patAutomaton\" should be specified\n");
        RAISE Scan.BadFormat
      END;
      
      o.verbose := PP.KeywordPresent("-verbose");
      
      IF PP.KeywordPresent("-breakChar") THEN
        o.breakChar := Text.GetChar(PP.GetNext(), 0);
      ELSE
        o.breakChar := '|';
      END;

      o.ignoreChars := Char.Spaces;
      o.commentChars := SET OF CHAR {'#'};
      o.ignoreChars := Char.Spaces;

      o.validChars := Char.Graphics + SET OF CHAR{'\200'..'\377'};
      o.validChars := o.validChars - o.commentChars - o.ignoreChars;
      o.validChars := o.validChars - SET OF CHAR{o.breakChar};

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

      IF o.breakChar IN o.validChars THEN
        Wr.PutText(stderr, "GetOptions: -validChars should not include -breakChar.\n");
        RAISE Scan.BadFormat
      END;

      PP.UnparsedTail();
      PP.EndParsing()
    EXCEPT
    | ParseParams.Error => 
        Wr.PutText(stderr, "usage: \n");
        Wr.PutText(stderr, "  patsplit \\\n");
        Wr.PutText(stderr, "    [-patAutomaton file.dmp | -patterns file.pats]\\\n");
        Wr.PutText(stderr, "    [-verbose]\\\n");
        Wr.PutText(stderr, "    [-breakChar char] \\\n");
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

PROCEDURE ErrorCannotOpen(file: TEXT) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(stderr, "cannot open file \"" & file & "\"\n");
    Wr.Flush(stderr);
    Process.Exit(1)
  END ErrorCannotOpen;

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
END PatSplit.
