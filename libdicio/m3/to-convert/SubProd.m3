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

(* Last modified on Thu Mar 19 12:38:58 PST 1992 by stolfi                  *)

MODULE SubProd EXPORTS Main;

IMPORT Rd, Wr, Text, Fmt, FileRd, FileWr, OSError, ParseParams, Thread, Scan;
IMPORT Basics, StringPrinter, Reduced, Encoding, PlainEncoding;
FROM Basics IMPORT INT, NAT, POS, BOOL;
FROM Reduced IMPORT State, NullState;
FROM Stdio IMPORT stderr;

TYPE
  Options = RECORD
      loadFile: TEXT;     (* Name of file to load  *)
      prodFile: TEXT;    (* Name of file with maximal subproducs, or "" if none *)
      minPrefs: NAT;          (* Minimum number of prefixes in prefix set *)
      minSuffs: NAT;          (* Minimum number of suffixes in suffix set *)
    END;

VAR o: Options;

PROCEDURE Main() =

  VAR encoding: Encoding.T := PlainEncoding.New();

  VAR aut: Reduced.T;

  <* FATAL Basics.Abort, Rd.Failure, Wr.Failure, Thread.Alerted *>
  BEGIN
    o := GetOptions();

    (* Load the automaton: *)
    <* ASSERT NOT Text.Empty(o.loadFile) *>
    Wr.PutText(stderr, 
      "\n=== loading automaton from " & o.loadFile & " ===\n\n"
    );
    WITH rd = FileRd.Open(o.loadFile) DO
      aut := Reduced.Load(rd);
    END;

    (* Print report *)
    WITH ct = aut.Count(ARRAY OF State{aut.Root()}) DO
      Reduced.PrintCounts(stderr, ct)
    END;

    (* Print counts and labels if so requested: *)
    IF NOT Text.Empty(o.prodFile) THEN
      Wr.PutText(stderr, 
        "\n=== printing pairs of states with common suffixes to " 
        & o.prodFile & " ===\n\n"
      );
      WITH pr = FileWr.Open(o.prodFile) DO
        EnumSubProds(pr, aut, encoding)
      END;
    ELSE
      EnumSubProds(NIL, aut, encoding)
    END;

    Wr.Flush(stderr);
  END Main;

CONST MaxStateSetSize = 2;

PROCEDURE CountCommonSuffs(
    aut: Reduced.T;
    VAR s: ARRAY [0..MaxStateSetSize-1] OF Reduced.State; 
    ns: POS;
  ): NAT =
  (* Counts suffixes that are common to all states s[0], s[1], ...s[ns-1]. 
     Destroys the contents of s. *)

  PROCEDURE SortS() =
    (* Sorts the vector "s[]", eliminates repeated states: *)
    BEGIN
      (* Sort states: *)
      FOR j := 1 TO ns - 1 DO
        VAR k: INT := j;
            r: Reduced.State := s[j];
        BEGIN
          WHILE k > 0 AND r < s[k-1] DO
            s[k] := s[k-1]; DEC(k)
          END;
          s[k] := r
        END
      END;
      (* Eliminate repeated ones: *)
      VAR nn: POS := 1;
          sl: Reduced.State := s[0];
      BEGIN
        FOR i := 1 TO ns-1 DO 
          IF s[i] # sl THEN
            sl := s[i];
            s[nn] := sl;
            INC(nn)
          END
        END;
        ns := nn
      END;
    END SortS;

  VAR t: ARRAY [0..MaxStateSetSize-1] OF Reduced.State;
      let: Basics.Letter := LAST(Basics.Letter);

  PROCEDURE AdvanceS() = 
    (* Find a common letter "let" to all states in "s[]".
       Put the destinations in "t[]", rests of "s[]" in "s[]".
       If there is no common letter, sets "let := FIRST(Basics.Letter)".
       *)      

    VAR i: NAT := 0;
    BEGIN
      WHILE i < ns AND let > FIRST(Basics.Letter) DO
        WITH si = s[i], ti = t[i] DO
          IF aut.HasArcs(si) THEN
            WITH ai = aut.Last(si) DO
              IF ai.letter < let THEN
                (* "let" won't do; try again with "ai.letter" *)
                let := ai.letter;
                i := 0;
              ELSIF ai.letter = let THEN
                (* So far so good: *)
                ti := ai.dest;
                si := aut.Rest(si);
                INC(i)
              ELSE
                (* Ignore this arc, try next arc from same "si": *)
                si := aut.Rest(si)
              END
            END
          ELSE
            let := FIRST(Basics.Letter);
            RETURN
          END
        END
      END
    END AdvanceS;

  VAR ncs: NAT := 0;
  
  BEGIN
    (* Count suffixes of descendants: *)
    LOOP
      SortS();
      IF ns = 1 THEN 
        ncs := ncs + aut.NSuffs(s[0]);
        RETURN ncs
      END;
      AdvanceS();
      IF let = FIRST(Basics.Letter) THEN EXIT END;
      ncs := ncs + CountCommonSuffs(aut, t, ns);
    END;
    (* Count empty suffix, if all states are final: *)
    FOR i := 0 TO ns-1 DO
      IF NOT aut.Final(s[i]) THEN RETURN ncs END
    END;
    RETURN ncs + 1
  END CountCommonSuffs;

PROCEDURE EnumSubProds(
    <*UNUSED*> wr: Wr.T;              (* May be NIL *)
    aut: Reduced.T;
    <*UNUSED*> encoding: Encoding.T;
  ) =

  CONST MaxHist = 1000;
  VAR pairs_ncs_hist: ARRAY [0..MaxHist] OF LONGREAL; (* Num state pairs per ncs *)
      nStates: LONGREAL; (* Num of states of automaton (not of the DAG!). *)
      nHard: LONGREAL;   (* Num of state pairs whose ncs was actually computed *)
      nGood: LONGREAL;   (* Num of state pairs whose ncs was big enough *)
      VAR s: ARRAY [0..MaxStateSetSize-1] OF Reduced.State;    
      npi, npj: NAT;

  <* FATAL Wr.Failure, Thread.Alerted *> 
  BEGIN
  
    FOR k := 0 TO LAST(pairs_ncs_hist) DO 
      pairs_ncs_hist[k] := 0.0d0 
    END;
  
    FOR i := 1 TO aut.Root() DO
      IF aut.NPrefs(i) > 0 THEN 
        nStates := nStates + 1.0d0 
      END
    END;

    FOR i := 1 TO aut.Root() - 1 DO
      npi := aut.NPrefs(i);
      IF npi > 0 AND aut.NSuffs(i) >= o.minSuffs THEN
        FOR j := i+1 TO aut.Root() DO 
          npj := aut.NPrefs(j);
          IF npj > 0 AND aut.NSuffs(j) >= o.minSuffs THEN
            nHard := nHard + 1.0d0;
            s[0] := i;
            s[1] := j;
            WITH 
              ncs = CountCommonSuffs(aut, s, 2),
              ih = MIN(ncs, LAST(pairs_ncs_hist))
            DO
              IF ih >= o.minSuffs THEN
                pairs_ncs_hist[ih] := pairs_ncs_hist[ih] + 1.0d0;
                nGood := nGood + 1.0d0
              END;
            END
          END
        END;
      END;
    END;

    WITH 
      nPairs = nStates * (nStates - 1.0d0)/2.0d0
    DO
      Wr.PutText(stderr, "states =        " & 
        Fmt.Pad(Fmt.LongReal(nStates, 0, Fmt.Style.Fix), 12) & 
        "\n"
      );
      Wr.PutText(stderr, "total pairs =   " &
        Fmt.Pad(Fmt.LongReal(nPairs, 0, Fmt.Style.Fix), 12) & 
        "\n"
      );
      Wr.PutText(stderr, "checked pairs = " & 
        Fmt.Pad(Fmt.LongReal(nHard, 0, Fmt.Style.Fix), 12) & 
        "\n"
      );
      Wr.PutText(stderr, "good pairs =    " & 
        Fmt.Pad(Fmt.LongReal(nGood, 0, Fmt.Style.Fix), 12) & 
        "\n"
      );
      Wr.PutText(stderr, "\n");
      Wr.PutText(stderr, " suffs           pairs\n");
      Wr.PutText(stderr, "======   =============\n");
      FOR ih := o.minSuffs TO LAST(pairs_ncs_hist) DO
        IF pairs_ncs_hist[ih] > 0.0d0 THEN
          Wr.PutText(stderr, Fmt.Pad(Fmt.Int(ih), 6));
          IF ih = LAST(pairs_ncs_hist) THEN 
            Wr.PutChar(stderr, '+')
          ELSE
            Wr.PutChar(stderr, ' ')
          END;
          Wr.PutText(stderr, "  ");
          Wr.PutText(stderr, 
            Fmt.Pad(Fmt.LongReal(pairs_ncs_hist[ih], 0, Fmt.Style.Fix), 12)
          );
          Wr.PutText(stderr, "\n");
        END
      END;
    END;
  
  END EnumSubProds;

<*UNUSED*> 
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
  <* FATAL Scan.BadFormat *>
  BEGIN
    WITH pp = NEW(ParseParams.T).init(stderr) DO

      IF pp.keywordPresent("-loadFile") THEN
        o.loadFile := pp.getNext()
      ELSE
        o.loadFile := ""
      END;

      IF pp.keywordPresent("-prodFile") THEN
        o.prodFile := pp.getNext()
      ELSE
        o.prodFile := ""
      END;

      IF pp.keywordPresent("-minSuffs") THEN
        o.minSuffs := pp.getNextInt()
      ELSE
        o.minSuffs := 2
      END;

    pp.finish();
    RETURN o
  END GetOptions;

BEGIN
  Main()
END SubProd.

