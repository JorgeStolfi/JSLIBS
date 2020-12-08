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

MODULE ElemProd EXPORTS Main;

IMPORT Rd, Wr, Text, Fmt, FileRd, FileWr, OSError, ParseParams, Thread, Scan;
IMPORT Basics, Reduced, Encoding, PlainEncoding;
FROM Basics IMPORT NAT, POS;
FROM Reduced IMPORT State, String, NullState;
FROM Stdio IMPORT stderr;

TYPE
  Options = RECORD
      dirAutomaton: TEXT;    (* Name of direct automaton's dumpfile *)
      revAutomaton: TEXT;    (* Name of reverse automaton's dumpfile  *)
      prodFile: TEXT;        (* Name of out file for elem prods, or "" if none *)
      minPrefs: POS;         (* List an elem prod only if it has at least this many prefs *)
      minSuffs: POS;         (* List an elem prod only if it has at least this many suffs *)
      minBytesSaved: POS;    (* List an elem prod only if it saves this many bytes *)
    END;

VAR o: Options;

PROCEDURE Main() =

  VAR dir, rev: Reduced.T;
      encoding: Encoding.T := PlainEncoding.New();

  <* FATAL Rd.Failure, Wr.Failure, Thread.Alerted *>
  BEGIN
    o := GetOptions();

    (* Load the automata: *)
    dir := ReadAutomaton(o.dirAutomaton);
    rev := ReadAutomaton(o.revAutomaton);

    (* Cound (and possibly print) related state pairs: *)
    IF NOT Text.Empty(o.prodFile) THEN
      Wr.PutText(stderr, 
        "\n=== printing pairs of related states to " 
        & o.prodFile & " ===\n\n"
      );
      WITH 
        pr = FileWr.Open(o.prodFile)
      DO
        EnumElemProds(pr, encoding, dir, rev)
      END;
    ELSE
      EnumElemProds(NIL, encoding, dir, rev)
    END;

    Wr.Flush(stderr);
  END Main;

PROCEDURE ReadAutomaton(file: TEXT): Reduced.T =
  <* FATAL Rd.Failure, Wr.Failure, Thread.Alerted *>
  BEGIN
    <* ASSERT NOT Text.Empty(file) *>
    Wr.PutText(stderr, 
      "\n=== loading automaton from " & file & " ===\n\n"
    );
    WITH 
      rd = FileRd.Open(file),
      aut = Reduced.Load(rd),
      ct = aut.Count(ARRAY OF State{aut.Root()})
    DO
      Reduced.PrintCounts(stderr, ct);
      Wr.PutText(stderr, "\n");
      RETURN aut
    END;
  END ReadAutomaton;

PROCEDURE EnumElemProds(
    wr: Wr.T; (* May be NIL *)
    encoding: Encoding.T;
    dir, rev: Reduced.T;
  ) =
  
  CONST MaxHist = 1000;
  VAR pairs_nw_hist: ARRAY [0..MaxHist] OF LONGREAL; (* Good pairs per nps*npt *)
      pairs_bs_hist: ARRAY [0..MaxHist] OF LONGREAL; (* Good pairs per bytes saved *)
      pairs_mp_hist: ARRAY [0..MaxHist] OF LONGREAL; (* Good pairs per MIN(nps,npt) *)
      nDirStates: LONGREAL; (* Num of states of "dir" (not of the DAG!). *)
      nRevStates: LONGREAL; (* Num of states of "rev" (not of the DAG!). *)
      nPairs: LONGREAL;  (* Num of state pairs  *)
      nGood: LONGREAL;   (* Num of state pairs whose ncs was big enough *)
      nps, npt: NAT;
      nls, nlt: NAT;

  <* FATAL Wr.Failure, Thread.Alerted *> 
  BEGIN
  
    FOR k := 0 TO LAST(pairs_nw_hist) DO 
      pairs_nw_hist[k] := 0.0d0; 
      pairs_bs_hist[k] := 0.0d0; 
      pairs_mp_hist[k] := 0.0d0;
    END;

  
    WITH
      s0 = dir.Root(),
      t0 = rev.Root(),
      rsPred = NEW(REF ARRAY OF Reduced.State, s0 + 1),
      rlPred = NEW(REF ARRAY OF Basics.Letter, s0 + 1),
      sPred = rsPred^,  (* sPred[i] is some immediate Succ-inverse of state i *)
      lPred = rlPred^   (* lPred[i] is the letter that takes sPred[i] to i *)
    DO

      FOR s := 1 TO s0 DO
        IF dir.NPrefs(s) > 0 THEN 
          nDirStates := nDirStates + 1.0d0 
        END
      END;

      FOR t := 1 TO t0 DO
        IF rev.NPrefs(t) > 0 THEN 
          nRevStates := nRevStates + 1.0d0 
        END
      END;

      (* Initialize sPred, lPred: *)
      lPred[s0] := FIRST(Basics.Letter);
      sPred[s0] := NullState;
      FOR s := dir.Root() TO 1 BY -1 DO
        IF dir.NPrefs(s) > 0 AND dir.HasArcs(s) THEN
          VAR i: State := s;
          BEGIN
            WHILE dir.HasArcs(i) DO
              WITH
                ai = dir.Last(i),
                ri = dir.Rest(i)
              DO
                sPred[ai.dest] := s;
                lPred[ai.dest] := ai.letter;
                i := ri
              END
            END
          END;
        END
      END;
      
      (* Enumerate pairs and check them: *)
      FOR s := 1 TO s0 DO
        nps := dir.NPrefs(s);
        IF nps > 0 THEN
          nls := dir.NPrefLetters(s);
          FOR t := 1 TO t0 DO
            npt := rev.NPrefs(t);
            IF npt > 0 THEN
              nlt := rev.NPrefLetters(t);
              nPairs := nPairs + 1.0d0;
              VAR sx := s;
                  tx := t;
              BEGIN
                WHILE sx # s0 AND tx # NullState DO
                  tx := rev.Succ(tx, lPred[sx]);
                  sx := sPred[sx];
                  <* ASSERT sx # NullState *>
                END;
                IF tx # NullState 
                AND rev.Final(tx) THEN
                  nGood := nGood + 1.0d0;
                  WITH 
                    nw = nps * npt,
                    ih_nw = MIN(nw, LAST(pairs_nw_hist)),
                    lp = nps * nlt + nls * npt,     (* Letters in product *)
                    bp = lp + nps * npt,            (* Bytes in product *)
                    bf = (nps + nls) + (npt + nlt), (* Bytes in factors *)
                    bs = bp - bf,                   (* Bytes saved *)
                    ih_bs = MAX(0, MIN(bs, LAST(pairs_bs_hist))),
                    mp = MIN(nps, npt),
                    ih_mp = MIN(mp, LAST(pairs_mp_hist))
                  DO
                    pairs_nw_hist[ih_nw] := pairs_nw_hist[ih_nw] + 1.0d0;
                    pairs_bs_hist[ih_bs] := pairs_bs_hist[ih_bs] + 1.0d0;
                    pairs_mp_hist[ih_mp] := pairs_mp_hist[ih_mp] + 1.0d0;
                    IF wr # NIL
                    AND nps >= o.minPrefs
                    AND npt >= o.minSuffs
                    AND bs >= o.minBytesSaved
                    THEN
                      EnumPrefsProduct(wr, encoding, dir, s, rev, t)
                    END
                  END;
                END
              END
            END
          END
        END
      END;
      
    END;
    
    WITH 
      nPairs = nDirStates * nRevStates
    DO
      Wr.PutText(stderr, "dir states =    " & 
        Fmt.Pad(Fmt.LongReal(nDirStates, 0, Fmt.Style.Fix), 12) & 
        "\n"
      );
      Wr.PutText(stderr, "rev states =    " & 
        Fmt.Pad(Fmt.LongReal(nRevStates, 0, Fmt.Style.Fix), 12) & 
        "\n"
      );
      Wr.PutText(stderr, "total pairs =   " &
        Fmt.Pad(Fmt.LongReal(nPairs, 0, Fmt.Style.Fix), 12) & 
        "\n"
      );
      Wr.PutText(stderr, "good pairs =    " & 
        Fmt.Pad(Fmt.LongReal(nGood, 0, Fmt.Style.Fix), 12) & 
        "\n"
      );
      Wr.PutText(stderr, "\n");
      PrintHist(stderr, "words",       "pairs", pairs_nw_hist);
      PrintHist(stderr, "bytes saved", "pairs", pairs_bs_hist);
      PrintHist(stderr, "min prefs",   "pairs", pairs_mp_hist);
    END;
  
  END EnumElemProds;
  
PROCEDURE PrintHist(
    wr: Wr.T;
    iTitle: TEXT; (* Title for index column *)
    cTitle: TEXT; (* Title for count column *)
    READONLY hist: ARRAY OF LONGREAL;
  ) =
  <* FATAL Wr.Failure, Thread.Alerted *> 
  BEGIN
    Wr.PutText(wr, Fmt.Pad(iTitle, 14));
    Wr.PutText(wr, "   ");
    Wr.PutText(wr, Fmt.Pad(cTitle, 14));
    Wr.PutChar(wr, '\n');
    Wr.PutText(wr, "==============   ==============\n");
    FOR ih := 0 TO LAST(hist) DO
      IF hist[ih] > 0.0d0 THEN
        Wr.PutText(wr, Fmt.Pad(Fmt.Int(ih), 6));
        IF ih = LAST(hist) THEN 
          Wr.PutChar(wr, '+')
        ELSE
          Wr.PutChar(wr, ' ')
        END;
        Wr.PutText(wr, "  ");
        Wr.PutText(wr, Fmt.Pad(Fmt.LongReal(hist[ih], 0, Fmt.Style.Fix), 12));
        Wr.PutText(wr, "\n");
      END
    END;
    Wr.Flush(wr);
  END PrintHist;

PROCEDURE EnumPrefsProduct(
    wr: Wr.T;
    e: Encoding.T;
    dir: Reduced.T;
    s: State;
    rev: Reduced.T;
    t: State;
  ) =

  PROCEDURE EnumRevPrefs (* : PrefixAction *) (READONLY wDir: String) =

    PROCEDURE PrintWord (* : PrefixAction *) (READONLY wRev: String) =
      <* FATAL Wr.Failure, Thread.Alerted *> 
      BEGIN
        FOR i := LAST(wDir) TO 0 BY -1 DO
          e.PrintLetter(wr, wDir[i])
        END;
        Wr.PutChar(wr, ':');
        FOR i := 0 TO LAST(wRev) DO
          e.PrintLetter(wr, wRev[i])
        END;
        Wr.PutChar(wr, '\n');
      END PrintWord;

    <* FATAL Basics.Abort *>
    BEGIN
      rev.EnumPrefs(t, PrintWord)
    END EnumRevPrefs;

  <* FATAL Wr.Failure, Thread.Alerted, Basics.Abort *> 
  BEGIN
    Wr.PutText(wr, "Prefs(" & Fmt.Int(s) & ", dir)");
    Wr.PutText(wr, " x ");
    Wr.PutText(wr, "Prefs(" & Fmt.Int(t) & ", rev)");
    Wr.PutText(wr, " = ");
    Wr.PutText(wr, Fmt.Int(dir.NPrefs(s)));
    Wr.PutText(wr, " x ");
    Wr.PutText(wr, Fmt.Int(rev.NPrefs(t)));
    Wr.PutText(wr, " = ");
    Wr.PutText(wr, Fmt.Int(dir.NPrefs(s)*rev.NPrefs(t)));
    Wr.PutText(wr, " = \n");
    dir.EnumPrefs(s, EnumRevPrefs);
    Wr.Flush(wr);
  END EnumPrefsProduct;

PROCEDURE GetOptions(): Options =
  VAR o: Options;
  <* FATAL Scan.BadFormat *>
  BEGIN
    WITH pp = NEW(ParseParams.T).init(stderr) DO

      pp.getKeyword("-dirAutomaton");
      o.dirAutomaton := pp.getNext();

      pp.getKeyword("-revAutomaton");
      o.revAutomaton := pp.getNext();

      IF pp.keywordPresent("-prodFile") THEN
        o.prodFile := pp.getNext()
      ELSE
        o.prodFile := ""
      END;

      IF pp.keywordPresent("-minPrefs") THEN
        o.minPrefs := pp.getNextInt(1, LAST(POS))
      ELSE
        o.minPrefs := 1
      END;

      IF pp.keywordPresent("-minSuffs") THEN
        o.minSuffs := pp.getNextInt(1, LAST(POS))
      ELSE
        o.minSuffs := 1
      END;

      IF pp.keywordPresent("-minBytesSaved") THEN
        o.minBytesSaved := pp.getNextInt(1, LAST(POS))
      ELSE
        o.minBytesSaved := 1
      END;

    pp.finish();
    RETURN o
  END GetOptions;

BEGIN
  Main()
END ElemProd.

