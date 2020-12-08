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

MODULE TestReduced EXPORTS Main;

IMPORT Rd, Wr, Text, Fmt, FileRd, FileWr, OSError, ParseParams, Process, Thread;
IMPORT Basics, StringPrinter, Reduced, Encoding, PlainEncoding;
FROM Basics IMPORT NAT, BOOL;
FROM Basics IMPORT Done;
FROM Reduced IMPORT String, State, NullState;
FROM Stdio IMPORT stdin, stdout, stderr;

CONST
  Help = 
    "TestReduced \\\n" &
    "  [ -load <file> ] \\\n" &
    "  [ -dump <file> ] \\\n" &
    "  [ -words <file> ] \\\n" &
    "  [ -print <file> ] \\\n" &
    "  [ -counts <file> ] \\\n" &
    "  [ -effs <file> ] \\\n" &
    "  [ -prefsuff <file> ] \\\n";

TYPE
  Options = RECORD
      loadFileName: TEXT;     (* Name of file to load, or "" to build from stdin  *)
      dumpFileName: TEXT;     (* Name of dump file, or "" if none *)
      wordsFileName: TEXT;    (* Name of file for output wordlist, or "" if none  *)
      printFileName: TEXT;    (* Name of legible printout file, or "" if none *)
      countsFileName: TEXT;   (* Name of file for pref/suff/word cts, or "" if none *)
      effsFileName: TEXT;     (* Name of file for state efficiencies, or "" if none *)
      prefsuffFileName: TEXT; (* Name of file for state prefs&suffs, or "" if none *)
    END;

EXCEPTION
  MissingFinalNewline;

VAR o: Options;

PROCEDURE Main() =

  VAR rd: Rd.T := stdin;
      encoding: Encoding.T := PlainEncoding.New();

  PROCEDURE NextString(
      VAR (*IO*) s: REF String; 
      VAR (*OUT*) len: NAT;
      VAR (*OUT*) add: BOOL;
    ) RAISES {Done} =

    VAR c: CHAR;

    <* FATAL Rd.Failure, Thread.Alerted, MissingFinalNewline *>
    BEGIN
      LOOP
        TRY
          TRY 
            REPEAT c := Rd.GetChar(rd) UNTIL c # ' ';
          EXCEPT 
            Rd.EndOfFile => RAISE Done 
          END;
          IF c = '+' THEN
            add := TRUE
          ELSIF c = '-' THEN
            add := FALSE
          ELSE
            add := TRUE;
            Rd.UnGetChar(rd)
          END;
          encoding.ReadString(rd, s, len);
          RETURN
        EXCEPT
        | Encoding.BadChars(msg) =>
          Wr.PutText(stderr,
            "** Encoding.BadChars: " & msg & 
            " at byte " & Fmt.Int(Rd.Index(rd)) &
            ", line ignored\n"
          );
          TRY
            WHILE Rd.GetChar(rd) # '\n' DO (* Ok *) END
          EXCEPT
            Rd.EndOfFile => RAISE MissingFinalNewline
          END;
        END
      END
    END NextString;

  VAR aut: Reduced.T;

  <* FATAL Basics.Abort, OSError.E, Wr.Failure, Thread.Alerted *>
  BEGIN
    o := GetOptions();

    (* Load or initialize the automaton: *)
    IF Text.Empty(o.loadFileName) THEN
      Wr.PutText(stderr, "\n=== initializing the automaton ===\n\n");
      aut := Reduced.New(size := 90000);
      aut.doc := 
        "\n" &
        "  Nature and nature's laws lay hidden in night,   \n" &
        "  God said, `Let Newton Be!', and all was light.  \n" &
        "                           --Alexander Pope       \n" &
        "                                                  \n" &
        "  It did not last; the Devil howling, `Ho,        \n" &
        "  let Einstein be!', restored the status quo.     \n" &
        "                           --J. C. Squire         \n" &
        "\n";
    ELSE
      Wr.PutText(stderr, 
        "\n=== loading automaton from " & o.loadFileName & " ===\n\n"
      );
      WITH rd = FileRd.Open(o.loadFileName) DO
        aut := Reduced.Load(rd);
      END;
      (* Print report *)
      WITH ct = aut.Count(ARRAY OF State{aut.Root()}) DO
        Reduced.PrintCounts(stderr, ct)
      END
    END;

    (* Process additions/deletions from stdin: *)
    BEGIN
      Wr.PutText(stderr, 
        "\n=== adding/subtracting words from stdin ===\n\n"
      );
      aut.Build(
        next := NextString,
        wr := stderr,
        reportInterval := 1000
      );
      aut.doc := aut.doc & 
        "\n" &
        "processed by TestReduced.m3\n";
    END;

    (* Dump if so requested: *)
    IF NOT Text.Empty(o.dumpFileName) THEN
      Wr.PutText(stderr, "\n=== dumping automaton to " & o.dumpFileName & " ===\n\n");
      WITH wr = FileWr.Open(o.dumpFileName) DO
        Reduced.Dump(wr, aut)
      END;
    END;

    (* Print reconstituted wordlist if so requested: *)
    IF NOT Text.Empty(o.wordsFileName) THEN
      Wr.PutText(stderr, 
        "\n=== printing reconstituted word list to " & o.wordsFileName & " ===\n\n"
      );
      WITH pr = FileWr.Open(o.wordsFileName) DO
        PrintAllWords(pr, aut, encoding)
      END;
    END;

    (* Print legible automaton if so requested: *)
    IF NOT Text.Empty(o.printFileName) THEN
      Wr.PutText(stderr, 
        "\n=== printing automaton to " & o.printFileName & " ===\n\n"
      );
      WITH pr = FileWr.Open(o.printFileName) DO
        Reduced.Print(pr, aut, encoding)
      END;
    END;

    (* Print counts and labels if so requested: *)
    IF NOT Text.Empty(o.countsFileName) THEN
      Wr.PutText(stderr, 
        "\n=== printing counts and full state labels to " 
        & o.countsFileName & " ===\n\n"
      );
      WITH pr = FileWr.Open(o.countsFileName) DO
        PrintAllCountsAndLabels(pr, aut, encoding)
      END;
    END;

    (* Print state efficiencies if so requested: *)
    IF NOT Text.Empty(o.effsFileName) THEN
      Wr.PutText(stderr, 
        "\n=== printing state efficiencies to " & o.effsFileName & " ===\n\n"
      );
      WITH pr = FileWr.Open(o.effsFileName) DO
        PrintAllEffsAndSuffs(pr, aut, encoding, maxChars := 50)
      END;
    END;

    (* Print all suffix and prefix sets if so requested: *)
    IF NOT Text.Empty(o.prefsuffFileName) THEN
      Wr.PutText(stderr, 
        "\n=== printing suffixes/prefixes to " & o.prefsuffFileName & " ===\n\n"
      );
      WITH pr = FileWr.Open(o.prefsuffFileName) DO
        PrintAllPrefsAndSuffs(pr, aut, encoding, maxChars := 500)
      END;
    END;

    Wr.Flush(stderr);
    Wr.Flush(stdout);
  END Main;

PROCEDURE PrintAllCountsAndLabels(
    wr: Wr.T; 
    aut: Reduced.T;
    encoding: Encoding.T;
  ) =

  CONST
    StateDigits = 6;
    CountDigits = 6;
    Sep = ":";

  PROCEDURE PrintStateCountsAndLabel (* : Reduced.StateAction *) (
      <*UNUSED*> len: NAT;
      s: State;
      <*UNUSED*> final: BOOL;
    ) RAISES {} =
    <* FATAL Wr.Failure, Thread.Alerted *>
    BEGIN
      (* Print prefixes: *)
      Wr.PutText(wr, Fmt.Pad(Fmt.Int(s), StateDigits));
      Wr.PutText(wr, "  ");
      WITH
        np = aut.NPrefs(s),
        ns = aut.NSuffs(s),
        nw = np * ns
      DO
        Wr.PutText(wr, Fmt.Pad(Fmt.Int(np), CountDigits));
        Wr.PutChar(wr, ' ');
        Wr.PutText(wr, Fmt.Pad(Fmt.Int(ns), CountDigits));
        Wr.PutChar(wr, ' ');
        Wr.PutText(wr, Fmt.Pad(Fmt.Int(nw), CountDigits));
      END;
      Wr.PutText(wr, "  ");

      Wr.PutText(wr, aut.FullLabel(s, e := encoding, sep := Sep));
      Wr.PutChar(wr, '\n');
      Wr.Flush(wr)
    END PrintStateCountsAndLabel;

  <* FATAL Basics.Abort, Wr.Failure, Thread.Alerted *>
  BEGIN
    aut.EnumStates(ARRAY OF State{aut.Root()}, enter := PrintStateCountsAndLabel);
    Wr.Flush(wr);
  END PrintAllCountsAndLabels;

PROCEDURE PrintAllEffsAndSuffs(
    wr: Wr.T; 
    aut: Reduced.T; 
    encoding: Encoding.T;
    maxChars: NAT;
  ) =

  CONST
    StateDigits = 6;
    CountDigits = 6;
    EffDecimals = 3;
    EffChars = 3 + 1 + EffDecimals;

  VAR
    suffPrinter := StringPrinter.New(
      wr := wr,
      encoding := encoding,
      sep := ", ",
      empty := "()",
      maxChars := maxChars,
      rightMargin := LAST(NAT)
    );

  PROCEDURE PrintStateEffAndSuffs (* : Reduced.StateAction *) (
      <*UNUSED*> len: NAT;
      s: State;
      <*UNUSED*> final: BOOL;
    ) RAISES {} =
    <* FATAL Wr.Failure, Thread.Alerted *>
    BEGIN
      (* Print prefixes: *)
      Wr.PutText(wr, Fmt.Pad(Fmt.Int(s), StateDigits));
      Wr.PutText(wr, "  ");
      WITH
        np = aut.NPrefs(s),
        ns = aut.NSuffs(s),
        nw = np * ns,
        eff = FLOAT(nw)/FLOAT(np + ns -1)
      DO
        Wr.PutText(wr, Fmt.Pad(Fmt.Int(np), CountDigits));
        Wr.PutChar(wr, ' ');
        Wr.PutText(wr, Fmt.Pad(Fmt.Int(ns), CountDigits));
        Wr.PutChar(wr, ' ');
        Wr.PutText(wr, Fmt.Pad(Fmt.Int(nw), CountDigits));
        Wr.PutChar(wr, ' ');
        Wr.PutText(wr, Fmt.Pad(Fmt.Real(eff, Fmt.Style.Fix, EffDecimals), EffChars));
      END;
      Wr.PutText(wr, "  { ");
      aut.PrintSuffs(s, suffPrinter);
      Wr.PutText(wr, " }\n");
      Wr.Flush(wr)
    END PrintStateEffAndSuffs;

  <* FATAL Basics.Abort, Wr.Failure, Thread.Alerted *>
  BEGIN
    aut.EnumStates(ARRAY OF State{aut.Root()}, enter := PrintStateEffAndSuffs);
    Wr.Flush(wr);
  END PrintAllEffsAndSuffs;

PROCEDURE PrintAllWords(
    wr: Wr.T;
    aut: Reduced.T;
    encoding: Encoding.T;
  ) =
(*
  Prints all words accepted by the automaton; also checks
  the "Rank" method (direct and reverse).
  *)

  VAR root: State := aut.Root();
      dirRank: NAT := 0;
      revRank: NAT := aut.NSuffs(root);

  PROCEDURE PrintWord (* : Reduced.SuffixAction *) (
      READONLY w: String;
    ) RAISES {} =
    <* FATAL Wr.Failure, Thread.Alerted *>
    BEGIN
      (* Print word *)
      FOR i := 0 TO LAST(w) DO encoding.PrintLetter(wr, w[i]) END;
      Wr.PutChar(wr, '\n');
      (* Update rank counters, and check "Rank" method: *)
      WITH r = aut.Rank(root, w) DO <* ASSERT r = dirRank *> END;
      INC(dirRank);
      DEC(revRank);
      WITH r = aut.Rank(root, w, reverse := TRUE) DO <* ASSERT r = revRank *> END;
      Wr.Flush(wr);
    END PrintWord;

  <* FATAL Basics.Abort, Wr.Failure, Thread.Alerted *>
  BEGIN
    aut.EnumSuffs(root, action := PrintWord);
    Wr.Flush(wr);
  END PrintAllWords;

PROCEDURE PrintAllPrefsAndSuffs(
    wr: Wr.T;
    aut: Reduced.T;
    encoding: Encoding.T;
    maxChars: NAT := LAST(NAT);
    lineWidth: NAT := 75;
  ) =
(*
  Prints the suffixes and prefixes of every state (truncating
  each at "maxChars"), in the format

     12323  Pref =   3 { abac, babac, abaremotem }
            Suff =   2 { o, os }

  *)

  CONST
    StateDigits = 6;
    LabelLength = 9;
    PrefLabel = "  Pref = ";
    SuffLabel = "  Suff = ";
    CountDigits = 6;
    OpenBrace = "{ ";
    CloseBrace = " }";
    BraceWidth = 2;
    LineIndent = StateDigits + LabelLength + CountDigits + 1 + BraceWidth;

  VAR spr := StringPrinter.New(
        wr := wr,
        encoding := encoding,
        empty := "()",
        sep := ", ",
        maxChars := maxChars,
        rightMargin := lineWidth,
        leftMargin := LineIndent,
        initialColumn := LineIndent
      );

  PROCEDURE PrintStatePrefsAndSuffs (* : Reduced.StateAction *) (
      <*UNUSED*> len: NAT;
      s: State;
      <*UNUSED*> final: BOOL;
    ) RAISES {} =
    <* FATAL Wr.Failure, Thread.Alerted *>
    BEGIN
      <* ASSERT s # NullState *>

      (* Print prefixes: *)
      Wr.PutText(wr, Fmt.Pad(Fmt.Int(s), StateDigits));
      Wr.PutText(wr, PrefLabel);
      WITH np = aut.NPrefs(s) DO
        Wr.PutText(wr, Fmt.Pad(Fmt.Int(np), CountDigits));
        Wr.PutChar(wr, ' ');
        IF np = 0 THEN
          Wr.PutText(wr, "{ }")
        ELSE
          Wr.PutText(wr, OpenBrace);
          spr.Reset();
          aut.PrintPrefs(s, spr);
          Wr.PutText(wr, CloseBrace);
        END;
      END;
      Wr.PutChar(wr, '\n');

      (* Print suffixes: *)
      FOR i := 0 TO StateDigits-1 DO Wr.PutChar(wr, ' ') END;
      Wr.PutText(wr, SuffLabel);
      WITH ns = aut.NSuffs(s) DO
        Wr.PutText(wr, Fmt.Pad(Fmt.Int(ns), CountDigits));
        Wr.PutChar(wr, ' ');
        IF ns = 0 THEN
          Wr.PutText(wr, "{ }")
        ELSE
          Wr.PutText(wr, OpenBrace);
          spr.Reset();
          aut.PrintSuffs(s, spr);
          Wr.PutText(wr, CloseBrace);
        END;
      END;
      Wr.PutChar(wr, '\n');
      Wr.PutChar(wr, '\n');
      Wr.Flush(wr)
    END PrintStatePrefsAndSuffs;

  <* FATAL Basics.Abort, Wr.Failure, Thread.Alerted *>
  BEGIN
    <* ASSERT LabelLength = Text.Length(PrefLabel) *>
    <* ASSERT LabelLength = Text.Length(SuffLabel) *>
    <* ASSERT BraceWidth = Text.Length(OpenBrace) *>
    <* ASSERT BraceWidth = Text.Length(CloseBrace) *>
    aut.EnumStates(ARRAY OF State{aut.Root()}, enter := PrintStatePrefsAndSuffs);
    Wr.Flush(wr);
  END PrintAllPrefsAndSuffs;

PROCEDURE GetOptions(): Options =
  VAR o: Options;
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    WITH pp = NEW(ParseParams.T).init(stderr) DO
      TRY

        IF pp.keywordPresent("-load") THEN
          o.loadFileName := pp.getNext()
        ELSE
          o.loadFileName := ""
        END;

        IF pp.keywordPresent("-dump") THEN
          o.dumpFileName := pp.getNext()
        ELSE
          o.dumpFileName := ""
        END;

        IF pp.keywordPresent("-words") THEN
          o.wordsFileName := pp.getNext();
        ELSE
          o.wordsFileName := ""
        END;

        IF pp.keywordPresent("-print") THEN
          o.printFileName := pp.getNext()
        ELSE
          o.printFileName := ""
        END;

        IF pp.keywordPresent("-counts") THEN
          o.countsFileName := pp.getNext()
        ELSE
          o.countsFileName := ""
        END;

        IF pp.keywordPresent("-effs") THEN
          o.effsFileName := pp.getNext()
        ELSE
          o.effsFileName := ""
        END;

        IF pp.keywordPresent("-prefsuff") THEN
          o.prefsuffFileName := pp.getNext()
        ELSE
          o.prefsuffFileName := ""
        END;

        pp.finish();
      EXCEPT
      | ParseParams.Error =>
          Wr.PutText(stderr, "usage: " & Help);
          Process.Exit(1);
      END
    END;    
    RETURN o
  END GetOptions;

BEGIN
  Main()
END TestReduced.

