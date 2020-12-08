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

MODULE TestReducedPair EXPORTS Main;

IMPORT
  Wr, Thread, Text, Fmt, FileRd, FileWr, OSError, Process, ParseParams,
  Basics, StringPrinter, Encoding, PlainEncoding, Reduced, ReducedPair;
FROM Basics IMPORT NAT;
FROM Basics IMPORT WrNat;
FROM Reduced IMPORT State;
FROM ReducedPair IMPORT BoolOp, Which;
FROM Stdio IMPORT stdout, stderr;

CONST
  Help =
    "TestReducedPair \\\n" &
    "  -load <file> \\\n" &
    "  [ -bools <file> ] \\\n" &
    "  [ -paths <file> ] \\\n";

TYPE
  Options = RECORD
      load: TEXT;   (* Name of file to load  *)
      bools: TEXT;  (* Name of file for BoolOps & suffs, or "" if none *)
      paths: TEXT;  (* Name of file for paths in product, or "" if none *)
    END;

VAR o: Options;

PROCEDURE Main() =

  VAR aut: Reduced.T;
      encoding: Encoding.T := PlainEncoding.New();

  <* FATAL Wr.Failure, Thread.Alerted, OSError.E *>
  BEGIN
    o := GetOptions();

    (* Build or load the automaton: *)
    <* ASSERT NOT Text.Empty(o.load) *>
    BEGIN
      Wr.PutText(stderr, 
        "\n=== loading automaton from " & o.load & " ===\n\n"
      );
      WITH rd = FileRd.Open(o.load) DO
        aut := Reduced.Load(rd)
      END;
      (* Print report *)
      WITH ct = aut.Count(ARRAY OF State{aut.Root()}) DO
        Wr.PutText(stderr, "root state = " & Fmt.Pad(Fmt.Int(aut.Root()), 8) & "\n");
        Wr.PutText(stderr, "states =     " & Fmt.Pad(Fmt.Int(ct.states), 8) & "\n");
        Wr.PutText(stderr, "sub-states = " & Fmt.Pad(Fmt.Int(ct.substates), 8) & "\n");
        Wr.PutText(stderr, "finals =     " & Fmt.Pad(Fmt.Int(ct.finals), 8) & "\n");
        Wr.PutText(stderr, "arcs =       " & Fmt.Pad(Fmt.Int(ct.arcs), 8) & "\n");
      END;
    END;

    (* Print all bool ops of all suffix sets if so requested: *)
    IF NOT Text.Empty(o.bools) THEN
      Wr.PutText(stderr, 
        "\n=== printing bool ops of suffix sets to " & o.bools & " ===\n\n"
      );
      WITH pr = FileWr.Open(o.bools) DO
        PrintAllBoolOps(pr, aut, encoding, maxChars := 500)
      END;
    END;

    (* Print all paths in product automaton if so requested: *)
    IF NOT Text.Empty(o.paths) THEN
      Wr.PutText(stderr, 
        "\n=== printing product paths to " & o.paths & " ===\n\n"
      );
      WITH pr = FileWr.Open(o.paths) DO
        PrintAllPairPaths(pr, aut, encoding)
      END;
    END;

    Wr.Flush(stderr);
    Wr.Flush(stdout);
  END Main;

PROCEDURE PrintAllBoolOps(
    wr: Wr.T;
    aut: Reduced.T;
    encoding: Encoding.T;
    maxChars: NAT := LAST(NAT);
    lineWidth: NAT := 75;
  ) =
(*
  Prints all boolean operations of the suffix sets of every
  pair of states (truncating each at "maxLetters"), in the format

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
      encoding := encoding,
      empty := "()",
      sep := ", ",
      maxChars := maxChars,
      rightMargin := lineWidth,
      leftMargin := LineIndent,
      initialColumn := LineIndent
    );

    pair: ReducedPair.T := ReducedPair.New(aut, aut);

  PROCEDURE PrintStateBoolOps (s: ReducedPair.State) RAISES {} =

    VAR argN: ARRAY Which OF NAT;
        opN: ARRAY BoolOp OF NAT;

    <* FATAL Wr.Failure, Thread.Alerted, Basics.Abort *>
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

  <* FATAL Wr.Failure, Thread.Alerted *>
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
    encoding: Encoding.T;
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
    <* FATAL Wr.Failure, Thread.Alerted *>
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
      <* FATAL Wr.Failure, Thread.Alerted *>
      BEGIN
        FOR i := 1 TO n DO Wr.PutChar(wr, ' ') END;
      END Indent;

    PROCEDURE PathEnter (* : ReducedPair.StateAction *) (
        len: NAT;
        s: ReducedPair.State;
        final: ReducedPair.Bools;
      ) RAISES {} =
      <* FATAL Wr.Failure, Thread.Alerted *>
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
      <* FATAL Wr.Failure, Thread.Alerted *>
      BEGIN
        Indent(len*4 + 4);
        Wr.PutText(wr, "push ");
        WrNat(wr, i);
        Wr.PutChar(wr, ':');
        encoding.PrintLetter(wr, arc.letter);
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
      <* FATAL Wr.Failure, Thread.Alerted *>
      BEGIN
        Indent(len*4 + 4);
        Wr.PutText(wr, "pop  ");
        WrNat(wr, i);
        Wr.PutChar(wr, ':');
        encoding.PrintLetter(wr, arc.letter);
        Wr.PutText(wr, " <- ");
        WITH d = arc.dest DO PrintStatePair(d, pair.Final(d)) END;
        Wr.PutChar(wr, '\n')
      END PathPop;

    PROCEDURE PathExit (* : ReducedPair.StateAction *) (
        len: NAT;
        s: ReducedPair.State;
        final: ReducedPair.Bools;
      ) RAISES {} =
      <* FATAL Wr.Failure, Thread.Alerted *>
      BEGIN
        Indent(len*4 + 2);
        Wr.PutText(wr, "exit  ");
        PrintStatePair(s, final);
        Wr.PutChar(wr, '\n')
      END PathExit;

    <* FATAL Wr.Failure, Thread.Alerted, Basics.Abort *>
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

  <* FATAL Wr.Failure, Thread.Alerted *>
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

PROCEDURE GetOptions(): Options =
  VAR o: Options;
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    TRY
      WITH pp = NEW(ParseParams.T).init(stderr) DO
        pp.getKeyword("-load");
        o.load := pp.getNext();

        IF pp.keywordPresent("-bools") THEN
          o.bools := pp.getNext()
        ELSE
          o.bools := ""
        END;

        IF pp.keywordPresent("-paths") THEN
          o.paths := pp.getNext()
        ELSE
          o.paths := ""
        END;

        pp.finish();
      END
    EXCEPT
    | ParseParams.Error =>
        Wr.PutText(stderr, "Usage: " & Help);
        Process.Exit(1)
    END;
    RETURN o
  END GetOptions;

BEGIN
  Main()
END TestReducedPair.

