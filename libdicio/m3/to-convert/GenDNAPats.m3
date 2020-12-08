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

MODULE GenDNAPats EXPORTS Main;

IMPORT 
  Rd, Wr, Thread, Fmt, Text, Word, Scan, Process,
  RandomPerm, Random;
IMPORT ParseParams AS PP;
IMPORT Basics;
FROM Basics IMPORT NAT;
FROM Stdio IMPORT stdout, stderr;

TYPE 
  Options = RECORD
    order: NAT;                  (* Length of patterns to generate *)
    alpha: NAT;                  (* Number of letters *)
    letters: REF ARRAY OF CHAR;  (* Alphabet *)
    goodEdges: BOOLEAN;          (* TRUE outputs the good edges, FALSE for cut edges *)
    verbose: BOOLEAN;            (* TRUE to mumble while thinking *)
  END;
  
TYPE
  StringCode = NAT;
  (* 
    A string  $x=(x_0,\dd x_n)$is represented by its StringCode.
    which is its $m$-adic value; that is, the number
    $\rank{x} = \sum (i \mapsto x_i m^i)$.  Note that $x_0$ 
    is the {\em least} significant letter. 
  *)

PROCEDURE Main() =
  <* FATAL Rd.Failure, Wr.Failure, Thread.Alerted *>
  BEGIN
    WITH
      o = GetOptions(),
      coins = NEW(Random.Default).init(TRUE),
      erep  =  (* Edge replement table *)
        ReplementTable(o.alpha, o.order)^,
      erot  =  (* Table of cyclic LeftShift for edges *)
        RotTable(o.alpha, o.order)^,
      etor  =  (* Table of cyclic Rightshift for edges *)
        PermInverse(erot)^,
      vrep  =  (* Vertex replement table *)
        ReplementTable(o.alpha, o.order-1)^,

      ecut = (* Current edge cutset (one per basic cycle) *)
        BasicCycles(o.alpha, o.order, erot, erep)^
    DO
      IF o.verbose THEN 

        Wr.PutText(stderr, "Vertex replement table:\n");
        WriteStringTable(stderr, o.alpha, o.order-1, vrep,  o.letters^);

        Wr.PutText(stderr, "Edge replement table:\n");
        WriteStringTable(stderr, o.alpha, o.order,   erep,  o.letters^);

        Wr.PutText(stderr, "Edge left rotation (rot) table:\n");
        WriteStringTable(stderr, o.alpha, o.order,   erot,  o.letters^);

        Wr.PutText(stderr, "Edge right rotation (tor) table:\n");
        WriteStringTable(stderr, o.alpha, o.order,   etor,  o.letters^);


        Wr.PutText(stderr, "Initial cutset:\n");
        WriteStringList(stderr, o.alpha, o.order,   ecut,  o.letters^);

      END;

      WritePatterns(stdout, 
        o.alpha, o.order,   
        pats := ecut, rept := erep,
        letters := o.letters^, 
        cutMark := NOT o.goodEdges
      );

      IF o.order MOD 2 = 0 AND o.alpha MOD 2 = 0 THEN
        WITH
          vrank =  (* Rank of vertices *)
            PickSymmetricRank(o.alpha, o.order-1, vrep, coins)^,
          vknar =  (* Inverse of rank: ordered list of StringCodes *)
            PermInverse(vrank)^
        DO
          IF o.verbose THEN
            Wr.PutText(stderr, "Initial vertex ordering:\n");
            WriteStringList(stderr, o.alpha, o.order-1, vknar, o.letters^);
          END;
        END;
      END
    END;
    
    Wr.Close(stdout);
    Wr.Close(stderr);
  END Main;
  
PROCEDURE ReplementTable(
    m: NAT; (* Alphabet size *)
    n: NAT; (* String length *)
  ): REF ARRAY OF StringCode =
  (* 
    Returns a table that maps each StringCode in $\Sigma_m^n$ 
    to its replement (alphabet complement of string reversal).
  *)
  VAR nStrings: NAT := 1;
  BEGIN
    WITH
      r = NEW(REF ARRAY OF StringCode, Expt(m, n)),
      rept = r^
    DO
      rept[0] := 0;
      FOR k := 1 TO n DO 
        (* Compute the table for $k$ letters from that of $k-1$ letters: *)
        FOR x := 0 TO nStrings-1 DO 
          FOR s := m-1 TO 0 BY -1 DO
            rept[s * nStrings + x] := rept[x] * m + (m - 1 - s);
          END
        END;
        nStrings := m * nStrings
      END;
      RETURN r
    END
  END ReplementTable;

PROCEDURE RotTable(
    m: NAT; (* Alphabet size *)
    n: NAT; (* String length *)
  ): REF ARRAY OF StringCode =
  (* 
    Returns a table that maps each StringCode in $\Sigma_m^n$ 
    to its cyclic left shift (Rot).
  *)
  BEGIN
    WITH
      ns = Expt(m, n),
      ns1 = ns DIV m,
      r = NEW(REF ARRAY OF StringCode, ns),
      rot = r^
    DO
      FOR i := 0 TO ns - 1 DO 
        rot[i] := (i DIV m) + (i MOD m) * ns1
      END;
      RETURN r
    END
  END RotTable;

PROCEDURE PickSymmetricRank(
    m: NAT;  (* Alphabet size *)
    n: NAT;  (* String length *)
    READONLY rept: ARRAY OF StringCode;  (* Replement table *)
    coins: Random.T;            (* Source of randomness *)
  ): REF ARRAY OF NAT =
  (* 
    Returns a table of random ranks for the vertices, 
    symmetric under replement: 
  *)
  VAR k: NAT := 0;
  BEGIN
    <* ASSERT m MOD 2 = 0 *>
    <* ASSERT n MOD 2 = 1 *>
    WITH
      nv = Expt(m, n),
      rrank = NEW(REF ARRAY OF StringCode, nv),
      rank = rrank^,
      perm = RandomPerm.NewArr(coins, nv DIV 2)^
    DO
      FOR i := 0 TO LAST(rank) DO 
        WITH j = rept[i] DO 
          IF i < j THEN 
            IF coins.integer() MOD 2 = 0 THEN 
              rank[i] := perm[k]; rank[j] := LAST(rank) - perm[k];
            ELSE
              rank[j] := perm[k]; rank[i] := LAST(rank) - perm[k];
            END;
            INC(k)
          END
        END
      END;
      <* ASSERT k = nv DIV 2 *>
      RETURN rrank
    END;
  END PickSymmetricRank;

PROCEDURE BasicCycles(
    m: NAT; (* Alphabet Size *)
    n: NAT; (* String length *)
    READONLY cyct: ARRAY OF StringCode; (* A permutation of $\Sigma_m^n$ *)
    READONLY rept: ARRAY OF StringCode; (* An involution on $\Sigma_m^n$ *)
  ): REF ARRAY OF StringCode =
  (*
    Returns a list with one string of $\Sigma_m^n$ 
    from each basic cycle of the DeBruijn graph $G_m^n$.
    The $rept$ table is used only to select the representative
    string from each cycle.
  *)
  VAR k: NAT := 0;
  BEGIN
    WITH
      ns = Expt(m, n),
      nc = NumBasicCycles(m, n),
      r = NEW(REF ARRAY OF StringCode, nc),
      cycle = r^
    DO
      FOR i := 0 TO ns - 1 DO 
        WITH
          j = CycleRep(i, cyct, rept) 
        DO
          IF i = j THEN 
            cycle[k] := i; INC(k)
          END
        END
      END;
      <* ASSERT k = nc *>
      RETURN r
    END
  END BasicCycles;
  
PROCEDURE NumBasicCycles(m, n: NAT): NAT =
  (*
    Computes the number of basic cycles of $G_m^n$
    These are the orbits of $Rot$ in $\Sigma_m^n$.
    These are also the distinct periodic sequences whose period 
    divides $n$, including those with periods that divide $n$. 
  *)
  <* FATAL Wr.Failure, Thread.Alerted *>
  VAR e: NAT := 1;
  BEGIN
    WITH
      epc = NEW(REF ARRAY OF NAT, n+1)^, (* Num edges in all primitive cycles of $G_m^i$ *)
      nbc = NEW(REF ARRAY OF NAT, n+1)^  (* Num basic cycles of $G_m^i$ *)
    DO
      FOR i := 1 TO n DO e := e*m; epc[i] := e; nbc[i] := 0 END;
      FOR i := 1 TO n DO 
        <* ASSERT epc[i] MOD 1 = 0 *>
        WITH npc = epc[i] DIV i DO
          (* npc is the number of cycles of length $i$ in $G_m^i$. *)
          nbc[i] := nbc[i] + npc;
          Wr.PutText(stderr, "m = " & Fmt.Int(m) & "  n = " & Fmt.Int(i));
          Wr.PutText(stderr, "  prim cycles = " & Fmt.Int(npc));
          Wr.PutText(stderr, "  basic cycles = " & Fmt.Int(nbc[i]));
          Wr.PutText(stderr, "\n");
          IF i = n THEN RETURN nbc[i] END;
          FOR j := 2*i TO n BY i DO 
            epc[j] := epc[j] - npc*i;
            nbc[j] := nbc[j] + npc;
          END;
        END
      END
    END 
  END NumBasicCycles;

PROCEDURE CycleRep(
    s: StringCode;
    READONLY cyct: ARRAY OF StringCode; (* Some permutation of $\Sigma_m^n$ *)
    READONLY rept: ARRAY OF StringCode; (* Some involution on $\Sigma_m^n$ *)
  ): StringCode =
  (*
    Returns a canonical representative $r(s)$ of the orbit of $s$ under $cyct$.
    Assuming that $rept$ is an involution that 
    takes cycles to cycles, the result satisfies $rept[r(s)] = r(rept[s])$,
    unless $s$ and $rept[s]$ belong to the same orbit, and that orbit contains
    no fixpoints of $rept$.
  *)
  VAR ruMin, rvMin: StringCode;
      uMin: NAT := NUMBER(rept);
      vMin: NAT := 2*NUMBER(rept);
      r: StringCode := s;
  BEGIN
    REPEAT
      IF rept[r] = r THEN
        WITH u = r DO
          IF u < uMin THEN ruMin := r; uMin := u END
        END
      ELSE
        WITH v = MIN(2 * r, 2 * rept[r] + 1) DO
          IF v < vMin THEN rvMin := r; vMin := v END
        END
      END;
      r := cyct[r]
    UNTIL r = s;
    IF uMin < NUMBER(rept) THEN RETURN ruMin ELSE RETURN rvMin END;
  END CycleRep;

PROCEDURE PermInverse(READONLY perm: ARRAY OF NAT): REF ARRAY OF NAT =
  BEGIN
    WITH
      r = NEW(REF ARRAY OF NAT, NUMBER(perm)),
      mrep = r^
    DO
      FOR i := 0 TO LAST(perm) DO 
        mrep[perm[i]] := i;
      END;
      RETURN r
    END;
  END PermInverse;

PROCEDURE WriteStringTable(
    wr: Wr.T;
    m: NAT;                               (* Alphabet size *)
    n: NAT;                               (* String length *)
    READONLY table: ARRAY OF StringCode;  (* Table to print *)
    READONLY letters: ARRAY OF CHAR;      (* Alphabet symbols *)
  ) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    FOR i := 0 TO LAST(table) DO 
      WriteString(wr, i, m, n, letters);
      Wr.PutText(wr, " -> ");
      WriteString(wr, table[i], m, n, letters);
      Wr.PutText(wr, "\n");
      Wr.Flush(wr);
    END;
  END WriteStringTable;

PROCEDURE WriteStringList(
    wr: Wr.T;
    m: NAT;                               (* Alphabet size *)
    n: NAT;                               (* String length *)
    READONLY table: ARRAY OF StringCode;  (* Table to print *)
    READONLY letters: ARRAY OF CHAR;      (* Alphabet symbols *)
  ) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    FOR i := 0 TO LAST(table) DO 
      Wr.PutText(wr, "[");
      Wr.PutText(wr, Fmt.Pad(Fmt.Int(i), 6));
      Wr.PutText(wr, "] = ");
      WriteString(wr, table[i], m, n, letters);
      Wr.PutText(wr, "\n");
      Wr.Flush(wr);
    END;
  END WriteStringList;

PROCEDURE WritePatterns(
    wr: Wr.T;
    m: NAT;                               (* Alphabet size *)
    n: NAT;                               (* String length *)
    READONLY pats: ARRAY OF StringCode;   (* Pattern strings *)
    READONLY rept: ARRAY OF StringCode;   (* String replement table *)
    READONLY letters: ARRAY OF CHAR;      (* Alphabet symbols *)
    cutMark: BOOLEAN;                     (* TRUE to insert cut mark *)
  ) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  VAR cutPos: NAT;
  BEGIN
    FOR i := 0 TO LAST(pats) DO 
      WITH p = pats[i] DO
        IF cutMark THEN 
          cutPos := n DIV 2;
          IF n MOD 2 = 1 AND rept[p] < p THEN 
            cutPos := n - cutPos
          END;
        ELSE
          cutPos := LAST(NAT)
        END;
        WritePattern(wr, p, m, n, letters, cutPos)
      END;
    END;
    Wr.PutText(stderr, Fmt.Pad(Fmt.Int(NUMBER(pats)), 8));
    Wr.PutText(stderr, " patterns written\n");
    Wr.Flush(stderr);
  END WritePatterns;

PROCEDURE WriteString(
    wr: Wr.T; 
    s: StringCode;
    m: NAT;
    n: NAT;
    READONLY letters: ARRAY OF CHAR;
  ) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    FOR i := 0 TO n-1 DO 
      WITH c = s MOD m DO
        Wr.PutChar(wr, letters[c]);
        s := s DIV m;
      END;
    END;
  END WriteString;

PROCEDURE WritePattern(
    wr: Wr.T;
    s: StringCode;                    (* Pattern string *)
    m: NAT;                           (* Alphabet size *)
    n: NAT;                           (* String length *)
    READONLY letters: ARRAY OF CHAR;  (* Alphabet symbols *)
    cutPos: NAT;                      (* Where to insert cut mark *)
  ) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    FOR i := 0 TO n-1 DO 
      IF i = cutPos THEN Wr.PutChar(wr, '|') END;
      WITH c = s MOD m DO
        Wr.PutChar(wr, letters[c]);
        s := s DIV m;
      END;
    END;
    IF n = cutPos THEN Wr.PutChar(wr, '|') END;
    Wr.PutChar(stdout, '\n');
    Wr.Flush(stdout);
  END WritePattern;

PROCEDURE Expt (b, e: NAT): NAT =
  VAR r: NAT := 1;
  BEGIN
    WHILE e > 0 DO r := r * b; DEC(e) END;
    RETURN r
  END Expt;
  
PROCEDURE GetOptions(): Options =

  VAR o: Options;

  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    TRY
      PP.BeginParsing(stderr);

      PP.GetKeyword("-order");
      o.order := PP.GetNextInt(2, Word.Size DIV 2);
      
      VAR lets: TEXT;
      BEGIN
        IF PP.KeywordPresent("-letters") THEN 
          lets := PP.GetNext()
        ELSE 
          lets := "ACGT";
        END;
        o.alpha := Text.Length(lets);
        o.letters := NEW(REF ARRAY OF CHAR, o.alpha);
        Text.SetChars(o.letters^, lets);
      END;
      o.goodEdges := PP.KeywordPresent("-goodEdges");
      o.verbose := PP.KeywordPresent("-verbose");

      PP.UnparsedTail();
      PP.EndParsing()
    EXCEPT
    | ParseParams.Error => 
        Wr.PutText(stderr, "usage: gendnapats -order nn \\\n");
        Wr.PutText(stderr, "    [-letters xxx] [-goodEdges]\n");
        Process.Exit(1);
    END;
    RETURN o
  END GetOptions;
  
BEGIN
  Main();
END GenDNAPats.

(*
      rank =  (* Rank of vertex in ordering *)
        ComputeVertexRank(o.order-1, o.alpha, rep, coins)^   
      WritePatterns(o.alpha, rank, o.letters);


PROCEDURE ComputeVertexRank(
    coins: Random.T; 
    READONLY rep: ARRAY OF NAT;
  ): REF ARRAY OF NAT =
  BEGIN
    WITH
      nv = NUMBER(rep),
      rank = NEW(REF ARRAY OF Word.T, nv),
      r = rank^
    DO
      FOR trial := 1 TO o.nTrials DO
        PickSymmetricRank(coins, rep, rank);
        VAR x: NAT := 0;
        BEGIN
          <* ASSERT nv DIV 2 = 0 *>
          WHILE x < nv DIV 2 DO
            
          END;
        END~

        (* Anneal the order so as to minimize the number of decreasing edges: *)
        FOR t := 0 TO o.nTries DO 
          WITH
            i = Word.And(coins.integer(), vmask)
          DO
            FOR c := 0 TO 3 DO 

            END~
          END~
        END~
      END;
      
      ~
    END~
      
        
    
  END ComputeVertexRank;

PROCEDURE Quality ()~
  
CONST 
  DNABases = ARRAY [0..3] OF CHAR{'A', 'C', 'G', 'T'};
  (* Complementary bases should have complementary indices. *)
  
PROCEDURE DNAComplement(w: Word.T; len: NAT): Word.T =
  (* Returns the complement of the base-reversal of the last "len" bases of "w". *)
  VAR r: Word.T := 0;
  BEGIN
    w := Word.Not(w);
    FOR i := 0 TO len-1 DO
      r := Word.Or(Word.Shift(r, 2), Word.And(w, 3));
      w := Word.Shift(w, -2)
    END;
    RETURN r
  END DNAComplement;

*)

