MODULE Huffes EXPORT Main;

IMPORT
  Rd, Wr, Char, Text, Fmt, FileStream, ParseParams,
  Basics, StringPrinter, Reduced, ReducedPair;
FROM Basics IMPORT INT, NAT, POS, BOOL;
FROM Basics IMPORT Done, Skip, Abort, WrNat;
FROM Reduced IMPORT Letter, String, State, Arc, NullState;
FROM ReducedPair IMPORT BoolOp, Which;
FROM Stdio IMPORT stdin, stdout, stderr;

TYPE
  Options = RECORD
      sep: CHAR;             (* Syllabe separator *)
      inputFileName: TEXT;   (* Name of input file (mandatory) *)
      freqsFileName: TEXT;   (* Name of output freq table file, or "" if none *)
    END;

VAR o: Options;

EXCEPTION 
  MissingFinalNewline;
  BadChar(CHAR); 
  BadLetter(Basics.Letter);
  EOF;  (* End-of-file; no more words *)
  EOL;  (* End-of-line; no more syllabes *)

PROCEDURE Main() =
  VAR code: Code.T;  (* A reduced automaton recognizing the code function. *)
  BEGIN
    o := GetOptions();

    (* Get the coding automaton: *)
    WITH rd = FileStream.OpenRead(o.codeFile) DO
      code := Reduced.Load(rd)
    END;
    
    (* Encode file: *)
    
    
    (* Print report *)
    WITH ct = aut.Count(ARRAY OF State{aut.Root()}) DO
      Wr.PutText(stderr, "root state = " & Fmt.Pad(Fmt.Int(aut.Root()), 8) & "\n");
      Wr.PutText(stderr, "states =     " & Fmt.Pad(Fmt.Int(ct.states), 8) & "\n");
      Wr.PutText(stderr, "sub-states = " & Fmt.Pad(Fmt.Int(ct.substates), 8) & "\n");
      Wr.PutText(stderr, "finals =     " & Fmt.Pad(Fmt.Int(ct.finals), 8) & "\n");
      Wr.PutText(stderr, "arcs =       " & Fmt.Pad(Fmt.Int(ct.arcs), 8) & "\n");
    END;

    Wr.Flush(stderr);
    Wr.Flush(stdout);
  END Main;
  
PROCEDURE CollectSyllabes(rd: Rd.T): Reduced.T RAISES {} =

  PROCEDURE NextSyllabe(VAR s: REF String; VAR len: NAT) RAISES {Done} =
    BEGIN
      ReadSyllabe(rd, s, len)
    END NextSyllabe;
    
  BEGIN
    WITH 
      aut = Reduced.Build(
        messages := stderr,
        next := NextSyllabe,
        estSize := 20000,
        reportInterval := 10000
      )
    DO
      RETURN aut
    END;
  END CollectSyllabes;

PROCEDURE ComputeFreqs(rd: Rd.T; syl: Reduced.T): REF ARRAY OF NAT RAISES {};
(*
  Returns a vector /freq/ where /freq[i]/ is th enumber of occurrences
  in /rd/ of the syllabe of rank /i/ in /syl/.
  Requires that the syllabes be all accepted by the /syl/.
  *)
  VAR s := NEW(REF String, 100);
      len: NAT;
  BEGIN
    WITH
      root = syl.Root(),
      freq = NEW(REF ARRAY OF NAT, syl.NSuffs(root)),
      f = freq^
    DO
      FOR i := 0 TO LAST(f) DO f[i] := 0 END;
      TRY
        LOOP
          REPEAT ReadSyllabe(rd, s, len) UNTIL len > 0;
          WITH 
            w = SUBARRAY(s^, 0, len),
            i = syl.Rank(root, w)
          DO
            IF NOT syl.Accepts(w) THEN RAISE InvaildSyllabe END;
            <* ASSERT i <= LAST(f) *>
            INC(f[i])
          END;
        END
      EXCEPT
      | Done => (* OK *)
      END;
      RETURN freq
    END;
  END ComputeFreqs;
  
PROCEDURE Encode(
    rd: Rd.T; 
    sep: CHAR; 
    syl: Reduced.T;
    code: Huffman.T;
    wr: Wr.T
  ) RAISES {} =
  VAR nWords: NAT := 0;
      sOut := NEW(String, 100);
      sSyl := NEW(String, 100);

  PROCEDURE ReadEncodeAndPrintWord(
    ) RAISES {EOF, BadChar} =
    VAR lenOut: NAT := 0;
        lenSyl: NAT;
    BEGIN
      TRY
        LOOP
          ReadSyllabe(rd, sep, (*IO*) sSyl, (*OUT*) lenSyl);
          WITH 
            rank = syl.Rank(syl.Root(), SUBARRAY(s^, 0, len))
          DO
            code.Append(rank, (*IO*) sOut, (*OUT*) lenOut)
          END;
        END
      EXCEPT
      | EOL => (* Ok *)
      END;
      (* Now print encoded string *)
      WITH w = SUBARRAY(sOut^, 0, lenOut) DO
        FOR i := 0 TO lenOut-1 DO
          PrintLetter(wr, w[i])
        END
      END;
      Wr.Flush(wr);
    END ReadEncodeAndPrintWord;

  BEGIN
    TRY
      LOOP
        TRY
          ReadEncodeAndPrintWord();
          INC(nWords)
        EXCEPT 
          BadChar(c) => ProcessBadChar(c, rd)
        END
      END
    EXCEPT
    | EOF => (* Ok *)
    END;
  END Encode;
  
PROCEDURE ReadSyllabe(
    rd: Rd.T;
    sep: CHAR;
    VAR (*IO*) s: REF String; 
    VAR (*OUT*) len: NAT
  ) RAISES {EOF, EOL, BadChar} =
(*
  Reads next syllabe from /rd/, converts each character
  to a Letter (using CharToLetter) and returns the result in /s, len/.
  A syllabe is a maximal sequence (possibly empty)
  of characters other than /sep/ or newline.  
  
  Raises EOL if the next character is '\n' (discards the character).
  Raises EOF if the reader is exhausted.
  *)
  VAR c: CHAR;
  BEGIN
    LOOP
      len := 0;
      TRY c := Rd.GetChar(rd) EXCEPT Rd.EndOfFile => RAISE EOF END;
      IF c = '\n' THEN RAISE EOL END;
      WHILE c # '\n' AND c # sep DO
        Basics.ExpandString(s, len+1);
        s[len] := CharToLetter(c);
        INC(len);
        TRY c := Rd.GetChar(rd) EXCEPT Rd.EndOfFile => RAISE MissingFinalNewline END
      END;
      IF len # 0 THEN RETURN END
    END
  END ReadSyllabe;

PROCEDURE ProcessBadChar(rd: Rd.T; c: CHAR) =
  BEGIN  
    Wr.PutText(stderr, 
      "** bad character (octal " &
      Fmt.Pad(Fmt.Int(ORD(c), base := 8), 3, '0') &
      ") at byte " & 
      Fmt.Int(Rd.Index(rd)) & 
      ", line ignored\n"
    );
    TRY 
      WHILE Rd.GetChar(rd) # '\n' DO (* Ok *) END
    EXCEPT 
      Rd.EndOfFile => RAISE MissingFinalNewline 
    END;
  END ProcessBadChar;  

PROCEDURE PrintAllCountsAndLabels(wr: Wr.T; aut: Reduced.T) =

  CONST 
    StateDigits = 6;
    CountDigits = 6;
    Sep = ":";
    
  PROCEDURE PrintStateCountsAndLabel (* : Reduced.StateAction *) (
      <*UNUSED*> len: NAT;
      s: State;
      <*UNUSED*> final: BOOL;
    ) RAISES {} =
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
      
      Wr.PutText(wr, aut.FullLabel(s, lpr := PrintLetter, sep := Sep));
      Wr.PutChar(wr, '\n');
      Wr.Flush(wr)
    END PrintStateCountsAndLabel;

  BEGIN
    aut.EnumStates(ARRAY OF State{aut.Root()}, enter := PrintStateCountsAndLabel);
    Wr.Flush(wr);
  END PrintAllCountsAndLabels;

PROCEDURE PrintAllEffsAndSuffs(wr: Wr.T; aut: Reduced.T; maxLetters: NAT) =

  CONST 
    StateDigits = 6;
    CountDigits = 6;
    EffDecimals = 3;
    EffChars = 3 + 1 + EffDecimals;
  
  VAR 
    suffPrinter := StringPrinter.New(
      wr := wr, 
      lpr := PrintLetter, 
      sep := ", ",
      empty := "()",
      maxLetters := maxLetters,
      maxLettersPerLine := LAST(NAT)
    );
  
  PROCEDURE PrintStateEffAndSuffs (* : Reduced.StateAction *) (
      <*UNUSED*> len: NAT;
      s: State;
      <*UNUSED*> final: BOOL;
    ) RAISES {} =
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
        Wr.PutText(wr, Fmt.Pad(Fmt.Real(eff, EffDecimals, Fmt.Style.Flo), EffChars));
      END;
      Wr.PutText(wr, "  { ");
      aut.PrintSuffs(s, suffPrinter);
      Wr.PutText(wr, " }\n");
      Wr.Flush(wr)
    END PrintStateEffAndSuffs;

  BEGIN
    aut.EnumStates(ARRAY OF State{aut.Root()}, enter := PrintStateEffAndSuffs);
    Wr.Flush(wr);
  END PrintAllEffsAndSuffs;

PROCEDURE PrintAllStrings(
    wr: Wr.T; 
    aut: Reduced.T;
    maxLetters: NAT := LAST(NAT);
    lineWidth: NAT := 75;
  ) =
(*
  Prints the suffixes and prefixes of every state (truncating 
  each at /maxLetters/), in the format
  
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
        lpr := PrintLetter, 
        empty := "()",
        sep := ", ",
        maxLetters := maxLetters,
        maxLettersPerLine := lineWidth - LineIndent,
        leftMargin := LineIndent
      );

  PROCEDURE PrintStateStrings (* : Reduced.StateAction *) (
      <*UNUSED*> len: NAT;
      s: State;
      <*UNUSED*> final: BOOL;
    ) RAISES {} =
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
    END PrintStateStrings;

  BEGIN
    <* ASSERT LabelLength = Text.Length(PrefLabel) *>
    <* ASSERT LabelLength = Text.Length(SuffLabel) *>
    <* ASSERT BraceWidth = Text.Length(OpenBrace) *>
    <* ASSERT BraceWidth = Text.Length(CloseBrace) *>
    aut.EnumStates(ARRAY OF State{aut.Root()}, enter := PrintStateStrings);
    Wr.Flush(wr);
  END PrintAllStrings;

PROCEDURE PrintAllBoolOps(
    wr: Wr.T; 
    aut: Reduced.T;
    maxLetters: NAT := LAST(NAT);
    lineWidth: NAT := 75;
  ) =
(*
  Prints all boolean operations of the suffix sets of every 
  pair of states (truncating each at /maxLetters/), in the format
  
     12323 op 3443
       
       Suff0 =   2 { abac, abaremotem }
       Suff1 =   2 { abac, babac }
       
       Union =   3 { abac, babac, abaremotem }
       Inter =   1 { abac }
       Diff  =   1 { babac }
       Symm  =   2 { babac, abaremotem }

  *)

  CONST 
    ArgLabel = ARRAY Which OF TEXT{
      "  Suff0 = ", 
      "  Suff1 = "
    };
    OpLabel = ARRAY BoolOp OF TEXT{
      "  Union = ", 
      "  Inter = ", 
      "  Diff  = ", 
      "  Symm  = "
    };
    LabelLength = 10;
    CountDigits = 6;
    OpenBrace = "{ ";
    CloseBrace = " }";
    BraceWidth = 2;
    LineIndent = LabelLength + CountDigits + 1 + BraceWidth;

  VAR 
    spr := StringPrinter.New(
      wr := wr, 
      lpr := PrintLetter, 
      empty := "()",
      sep := ", ",
      maxLetters := maxLetters,
      maxLettersPerLine := lineWidth - LineIndent,
      leftMargin := LineIndent
    );

    pair: ReducedPair.T := ReducedPair.New(aut, aut);
    
  PROCEDURE PrintStateBoolOps (s: ReducedPair.State) RAISES {} =

    VAR argN: ARRAY Which OF NAT;
        opN: ARRAY BoolOp OF NAT;
    
    BEGIN

      (* Print state *)
      Wr.PutText(wr, Fmt.Int(s[0]));
      Wr.PutText(wr, " op ");
      Wr.PutText(wr, Fmt.Int(s[1]));
      Wr.PutChar(wr, '\n');

      Wr.PutChar(wr, '\n');
      
      (* Print operands *)
      FOR i := 0 TO 1 DO 
        Wr.PutText(wr, ArgLabel[i]);
        WITH ns = aut.NSuffs(s[i]) DO
          argN[i] := ns;
          Wr.PutText(wr, Fmt.Pad(Fmt.Int(ns), CountDigits));
          Wr.PutChar(wr, ' ');
          IF ns = 0 THEN
            Wr.PutText(wr, "{ }")
          ELSE
            Wr.PutText(wr, OpenBrace);
            spr.Reset();
            aut.PrintSuffs(s[i], spr);
            Wr.PutText(wr, CloseBrace);
          END;
        END;
        Wr.PutChar(wr, '\n');
      END;

      Wr.PutChar(wr, '\n');

      (* Print operations: *)
      FOR op := FIRST(BoolOp) TO LAST(BoolOp) DO 
        Wr.PutText(wr, OpLabel[op]);
        WITH ns = pair.NSuffs(s, op) DO
          opN[op] := ns;
          Wr.PutText(wr, Fmt.Pad(Fmt.Int(ns), CountDigits));
          Wr.PutChar(wr, ' ');
          IF ns = 0 THEN
            Wr.PutText(wr, "{ }")
          ELSE
            Wr.PutText(wr, OpenBrace);
            spr.Reset();
            pair.PrintSuffs(s, op, spr);
            Wr.PutText(wr, CloseBrace);
          END;
        END;
        Wr.PutChar(wr, '\n');
      END;
      
      (* Check consistency of set sizes: *)
      <* ASSERT opN[BoolOp.Inter] <= MIN(argN[0], argN[1]) *>
      <* ASSERT opN[BoolOp.Union] = argN[0] + argN[1] - opN[BoolOp.Inter] *>
      <* ASSERT opN[BoolOp.Diff] = argN[0] - opN[BoolOp.Inter] *>
      <* ASSERT opN[BoolOp.Symm] = argN[0] + argN[1] - 2*opN[BoolOp.Inter] *>

      Wr.PutChar(wr, '\n');
      Wr.Flush(wr)
    END PrintStateBoolOps;

  BEGIN
    FOR i := 0 TO 1 DO 
      <* ASSERT LabelLength = Text.Length(ArgLabel[i]) *>
    END;
    FOR op := FIRST(BoolOp) TO LAST(BoolOp) DO 
      <* ASSERT LabelLength = Text.Length(OpLabel[op]) *>
    END;
    <* ASSERT BraceWidth = Text.Length(OpenBrace) *>
    <* ASSERT BraceWidth = Text.Length(CloseBrace) *>

    FOR s := 1 TO aut.Root() DO 
      FOR t := s+1 TO aut.Root() DO 
        IF aut.NPrefs(s) # 0 AND aut.NPrefs(t) # 0 THEN
          PrintStateBoolOps(ReducedPair.State{s, t})
        END
      END;
    END;

    Wr.Flush(wr);
  END PrintAllBoolOps;

PROCEDURE PrintAllPairPaths(
    wr: Wr.T; 
    aut: Reduced.T;
  ) =
(*
  Prints all proper paths in the product automaton, starting with
  every pair of states, as
  
     === pair ( 12323*, 3443 ) ===
     
       enter ( 12323*, 3343 )
         push 1:a -> ( 133, 145* )
         pop  1:a <- ( 133, 145* )
         push 2:c -> ( 999, 145* ) 
           enter ( 999, 145* )
           exit  ( 999, 145* )
         pop  2:c <- ( 999, 145* )
       exit ( 12323*, 3343 )
       
  *)
  
  PROCEDURE PrintStatePair(s: ReducedPair.State; f: ReducedPair.Bools) =
    BEGIN
      Wr.PutChar(wr, '(');
      Wr.PutText(wr, Fmt.Int(s[0]));
      IF f[0] THEN Wr.PutChar(wr, '*') END;
      Wr.PutChar(wr, ',');
      Wr.PutChar(wr, ' ');
      Wr.PutText(wr, Fmt.Int(s[1]));
      IF f[1] THEN Wr.PutChar(wr, '*') END;
      Wr.PutChar(wr, ')');
    END PrintStatePair;

  VAR
    pair: ReducedPair.T := ReducedPair.New(aut, aut);
    
  PROCEDURE PrintStatePairPaths (s: ReducedPair.State) RAISES {} =

    PROCEDURE Indent(n: NAT) =
      BEGIN
        FOR i := 1 TO n DO Wr.PutChar(wr, ' ') END;
      END Indent;

    PROCEDURE PathEnter (* : ReducedPair.StateAction *) (
        len: NAT; 
        s: ReducedPair.State;
        final: ReducedPair.Bools;
      ) RAISES {Skip, Abort} =
      BEGIN
        Indent(len*4 + 2);
        Wr.PutText(wr, "enter ");
        PrintStatePair(s, final);
        Wr.PutChar(wr, '\n')
      END PathEnter;

    PROCEDURE PathPush (* : ReducedPair.PathAction *) (
        len: NAT;
        <*UNUSED*> org: ReducedPair.State; 
        i: NAT; 
        arc: ReducedPair.Arc;
      ) RAISES {} =
      BEGIN
        Indent(len*4 + 4);
        Wr.PutText(wr, "push ");
        WrNat(wr, i);
        Wr.PutChar(wr, ':');
        PrintLetter(wr, arc.letter);
        Wr.PutText(wr, " -> ");
        WITH d = arc.dest DO PrintStatePair(d, pair.Final(d)) END;
        Wr.PutChar(wr, '\n')
      END PathPush;

    PROCEDURE PathPop (* : ReducedPair.PathAction *) (
        len: NAT;
        <*UNUSED*> org: ReducedPair.State; 
        i: NAT; 
        arc: ReducedPair.Arc;
      ) RAISES {} =
      BEGIN
        Indent(len*4 + 4);
        Wr.PutText(wr, "pop  ");
        WrNat(wr, i);
        Wr.PutChar(wr, ':');
        PrintLetter(wr, arc.letter);
        Wr.PutText(wr, " <- ");
        WITH d = arc.dest DO PrintStatePair(d, pair.Final(d)) END;
        Wr.PutChar(wr, '\n')
      END PathPop;

    PROCEDURE PathExit (* : ReducedPair.StateAction *) (
        len: NAT; 
        s: ReducedPair.State;
        final: ReducedPair.Bools;
      ) RAISES {Skip, Abort} =
      BEGIN
        Indent(len*4 + 2);
        Wr.PutText(wr, "exit  ");
        PrintStatePair(s, final);
        Wr.PutChar(wr, '\n')
      END PathExit;

    BEGIN
      (* Print state *)
      Wr.PutText(wr, "=== ");
      PrintStatePair(s, pair.Final(s));
      Wr.PutText(wr, " ===\n");

      Wr.PutChar(wr, '\n');
      
      pair.EnumPaths(
        s, 
        enter := PathEnter, 
        push := PathPush, 
        pop := PathPop, 
        exit := PathExit
      );

      Wr.PutChar(wr, '\n');
      Wr.Flush(wr)
    END PrintStatePairPaths;

  BEGIN

    FOR s := 1 TO aut.Root() DO 
      FOR t := s+1 TO aut.Root() DO 
        IF aut.NPrefs(s) # 0 AND aut.NPrefs(t) # 0 THEN
          PrintStatePairPaths(ReducedPair.State{s, t})
        END
      END;
    END;

    Wr.Flush(wr);
  END PrintAllPairPaths;

PROCEDURE PrintLetter(wr: Wr.T; letter: Basics.Letter) =
  BEGIN
    TRY
      Wr.PutChar(wr, LetterToChar(letter))
    EXCEPT
    | BadLetter(bad) => 
        Wr.PutChar(wr, '(');
        WrNat(wr, ORD(bad));
        Wr.PutChar(wr, ')');
    END;
  END PrintLetter;

PROCEDURE CharToLetter(c: CHAR): Reduced.Letter RAISES {BadChar} =
  BEGIN
    IF c IN Char.Controls THEN RAISE BadChar(c) END;
    RETURN ORD(c)
  END CharToLetter;
  
PROCEDURE LetterToChar(letter: Basics.Letter): CHAR RAISES {BadLetter} =
  BEGIN
    IF letter = 0 THEN RAISE BadLetter(letter) END;
    WITH c = VAL(letter, CHAR) DO
      IF c IN Char.Controls THEN RAISE BadLetter(letter) END;
      RETURN c
    END
  END LetterToChar;

PROCEDURE GetOptions(): Options =
  VAR o: Options;
  BEGIN
    ParseParams.BeginParsing(stderr);
      
      IF ParseParams.KeywordPresent("-print") THEN 
        o.printFileName := ParseParams.GetNext()
      ELSE
        o.printFileName := ""
      END;
    
    ParseParams.EndParsing();
    RETURN o
  END GetOptions;



BEGIN
  
END Huffes.
