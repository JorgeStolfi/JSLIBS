MODULE PermDAG;

IMPORT Code, DAG, Huffman, Util;
IMPORT OldFileFmt, FPut, NPut, FGet, NGet;
IMPORT FileRd, FileWr, Fmt, Time, Date, Rd, Text, Wr, Thread, OSError;

FROM Basics IMPORT Full;

FROM ReadOnlyPermDAG IMPORT 
  ComprHeader, ComprTrailer, DefaultCommentText, ExtendedLetter,
  FinalBitChar, InexistentLetter, Letter, State;

FROM Stdio IMPORT stderr;

<* FATAL Thread.Alerted, OSError.E, Rd.Failure, Wr.Failure *>

CONST
  PermDAGFormat : CARDINAL =  4;
  BackupPeriod  : Time.T = 3600.0d0; (* 3600 sec = 1 hour *)

TYPE
  Arc                   = DAG.Arc;

REVEAL
  T = Public BRANDED OBJECT
      comment            : TEXT;
      format             : CARDINAL;
      textperm           : TEXT;
      charperm           : REF ARRAY OF CHAR;
      chartolettertable  : ARRAY CHAR OF ExtendedLetter;
      root               : REF ARRAY OF State;
      dag                : DAG.T;
    OVERRIDES
      Accepts            := Accepts;
      AddSub             := AddSub;
      Conta              := Conta;
      Comment            := Comment;
      Crunch             := Crunch;
      Dump               := Dump;
      DumpCompr          := DumpCompr;
      DumpFixedBin       := DumpFixedBin;
      ExpandState        := ExpandState;
      Fold               := Fold;
      IncRoots           := IncRoots;
      IsFinal            := IsFinal;
      MakeTransition     := MakeTransition;
      Prm2Red            := Prm2Red;
      Root               := Root;
      SortPrm            := SortPrm;
      Spell              := Spell;
      Statistics         := Statistics;
      UnFold             := UnFold
    END;

TYPE
  AddressClass           = {Nul, Seq, Jmp};

CONST
  PermDAGFileType        : TEXT = "PermDAG.Dump";
  PermDAGFileVersion     : TEXT = "92-08-24";
  ReducedFileType        : TEXT = "Reduced.Dump";
  ReducedFileVersion     : TEXT = "91-12-21";
  FormatPrefix           : TEXT = "Format";
  RootPrefix             : TEXT = "PermDag.Root LAST";
VAR
  PermPrefix             : TEXT := "Permutation (" &
    Text.FromChar(FinalBitChar) & " = Char.NUL)";

CONST
  CharSize               : CARDINAL = BITSIZE(CHAR);
VAR
  BooleanSize            : CARDINAL := Code.BooleanSize();
  WordSize               : CARDINAL := Code.WordSize();
  WordSizeSize           : CARDINAL := Code.WordSizeSize() - 1;(* excludes 0 *)

PROCEDURE Accepts(permdag: T; s: State; word: TEXT): BOOLEAN =
  BEGIN
    WITH len = Text.Length(word),
         w   = word & Text.FromChar(FinalBitChar),
         dag = permdag.dag,
         t   = permdag.chartolettertable
         DO
      FOR i := 0 TO len DO
        WITH c = t[Text.GetChar(w, i)] DO
          WHILE s # 0 AND dag.Last(s).rd # c DO
            s := dag.Rest(s)
          END
        END;
        IF s = 0 THEN RETURN FALSE END;
        s := dag.Last(s).dest
      END;
      RETURN s = 0
    END
  END Accepts;
  
PROCEDURE AddSub(permdag: T; add: BOOLEAN; word: TEXT; level : CARDINAL := 0)
    RAISES {Full} =
  VAR
    ind : CARDINAL;
    lst : CARDINAL := Text.Length(word);
    wrd : REF ARRAY OF Letter := NEW(REF ARRAY OF Letter, 1 + lst);
    
  PROCEDURE VisitAdd(VAR s: State) RAISES {Full} =
    VAR
      arc  : Arc := Arc{rd := 0, dest := 0};
      rest : State := 0;
    BEGIN
      IF ind > lst THEN <* ASSERT s = 0 *> RETURN END;
      WITH let = wrd[ind],
           dag = permdag.dag
           DO
        IF s = 0 THEN arc.rd := let
        ELSE arc  := dag.Last(s); rest := dag.Rest(s)
        END;
        IF arc.rd = let THEN INC(ind); VisitAdd(arc.dest)
        ELSE VisitAdd(rest)
        END;
        s := dag.Append(arc, rest)
      END;
    END VisitAdd;
    
  PROCEDURE VisitSub(VAR s: State) RAISES {Full} =
    VAR
      arc  : Arc := Arc{rd := 0, dest := 0};
      rest : State := 0;
    BEGIN
      IF s = 0 THEN RETURN END;
      <* ASSERT ind <= lst *>
      WITH let = wrd[ind],
           dag = permdag.dag
           DO
        arc := dag.Last(s); rest := dag.Rest(s);
        IF arc.rd = let THEN
          INC(ind);
          VisitSub(arc.dest);
          IF arc.dest = 0 THEN s := rest; RETURN END;
        ELSE
          VisitSub(rest)
        END;
        s := dag.Append(arc, rest)
      END
    END VisitSub;
    
  BEGIN
    FOR i := 0 TO lst - 1 DO
      WITH c = Text.GetChar(word, i) DO
        <* ASSERT c # FinalBitChar *>
        wrd[i] := CharToLetter(permdag, c)
      END
    END;
    wrd[lst] := CharToLetter(permdag, FinalBitChar);
    <* ASSERT level <= LAST(permdag.root^) *>
    LOOP
      TRY
        ind := 0;
        IF add THEN VisitAdd(permdag.root[level])
        ELSE VisitSub(permdag.root[level])
        END;
        EXIT
      EXCEPT
        Full => 
          permdag.dag.Expand(2 + 11 * permdag.dag.MaxAllocState() DIV 10);
          Crunch(permdag)
      END
    END
  END AddSub;
  
PROCEDURE BuildPermDagTables(PermDag: T) =
  BEGIN
    WITH t = PermDag.textperm DO
      PermDag.charperm := NEW(REF ARRAY OF CHAR, Text.Length(t));
      Text.SetChars(PermDag.charperm^, t);
      FOR i := 0 TO LAST(PermDag.charperm^) DO
        PermDag.chartolettertable[PermDag.charperm[i]] := i
      END
    END
  END BuildPermDagTables;

PROCEDURE CharToLetter(permdag: T; c: CHAR): Letter =
  BEGIN
    WITH e = permdag.chartolettertable[c] DO
      IF e = InexistentLetter THEN
        WITH textperm = permdag.textperm,
             charperm = permdag.charperm
             DO
          e := Text.Length(textperm);
          textperm := textperm & Text.FromChar(c);
          charperm := NEW(REF ARRAY OF CHAR, 1 + e);
          Text.SetChars(charperm^, textperm)
        END
      END;
      RETURN e
    END
  END CharToLetter;
  
PROCEDURE ClassifyAddress(s: State; a: State): AddressClass =
  BEGIN
    IF a = 0     THEN RETURN AddressClass.Nul END;
    IF a = s - 1 THEN RETURN AddressClass.Seq END;
                      RETURN AddressClass.Jmp (* 1 <= a <= s - 2 *)
  END ClassifyAddress;
  
PROCEDURE Comment(permdag: T; comment: TEXT) =
  BEGIN
    permdag.comment := comment
  END Comment;
  
PROCEDURE Conta(permdag: T; wr: Wr.T) =
  TYPE
    Linha   = {i,n};
    Coluna  = {seq, pul, fim};
  CONST
    ltext   : ARRAY Linha OF TEXT = ARRAY Linha OF TEXT{"i", "n"};
    ctext   : ARRAY Coluna OF TEXT = ARRAY Coluna OF TEXT{"seq", "pul", "fim"};
  VAR
    nl      : CARDINAL := 0;
    nd      : CARDINAL := 0;
    linha   : Linha;
    coluna  : Coluna;
    matriz  : ARRAY Linha OF ARRAY Coluna OF CARDINAL := ARRAY Linha OF
      ARRAY Coluna OF CARDINAL{ARRAY Coluna OF CARDINAL{0, ..}, ..};
  BEGIN
    WITH perm    = permdag.charperm,
         dag     = permdag.dag,
         dagsize = dag.MaxState()
         DO
      FOR s := 1 TO dagsize DO
        WITH arc  = dag.Last(s), rest = dag.Rest(s) DO
          IF perm[arc.rd] = FinalBitChar THEN
            INC(nl);                       (* @ *)
            linha := Linha.i
          ELSE
            INC(nl); INC(nd);              (* r d *)
            IF arc.dest = s - 1 THEN
              linha := Linha.i
            ELSE
              linha := Linha.n
            END
          END;
          IF rest = 0 THEN
            INC(nl);                       (* | *)
            coluna := Coluna.fim
          ELSIF rest = s - 1 THEN
            coluna := Coluna.seq
          ELSE
            INC(nl); INC(nd);              (* > d *)
            coluna := Coluna.pul
          END
        END;
        INC(matriz[linha, coluna])
      END
    END;
    Wr.PutText(wr, "nl = " & Fmt.Int(nl) & " nd = " & Fmt.Int(nd) & 
      " total = " & Fmt.Int(nl + 2 * nd) & "\n");
    FOR l := FIRST(Linha) TO LAST(Linha) DO
      FOR c := FIRST(Coluna) TO LAST(Coluna) DO
        Wr.PutText(wr, ltext[l] & " " & ctext[c] & " " & 
          Fmt.Int(matriz[l,c]) & "\n")
      END
    END
  END Conta;

PROCEDURE Copy(frompermdag, topermdag: T) =
  BEGIN
    topermdag.format            := frompermdag.format;
    topermdag.comment           := frompermdag.comment;
    topermdag.textperm          := frompermdag.textperm;
    topermdag.charperm          := frompermdag.charperm;
    topermdag.chartolettertable := frompermdag.chartolettertable;
    topermdag.root              := frompermdag.root;
    topermdag.dag               := frompermdag.dag
  END Copy;
  
PROCEDURE Crunch(permdag: T) =
  BEGIN
    permdag.dag.Crunch(permdag.root^)
  END Crunch;
  
PROCEDURE Dump(permdag: T; wr: Wr.T) =
  BEGIN
    OldFileFmt.WriteHeader(wr, PermDAGFileType, PermDAGFileVersion);
    
    NPut.Int(wr, FormatPrefix, permdag.format);
    FPut.EOL(wr);
    
    OldFileFmt.WriteComment(wr, permdag.comment, '|');
    
    (* Careful - blanks are significant here! *)
    Wr.PutText(wr, PermPrefix);
    Wr.PutText(wr, " = ");
    Wr.PutText(wr, permdag.textperm);
    FPut.EOL(wr);
    
    WITH 
      root = permdag.root,
      Last = LAST(root^)
    DO
      NPut.Int(wr, RootPrefix, Last);
      FPut.Colon(wr);
      FPut.Int(wr, root[0]); 
      FOR i := 1 TO Last DO FPut.Space(wr); FPut.Int(wr, root[i]) END
    END;
    FPut.EOL(wr);
    
    DAG.Dump  (wr, permdag.dag);
    
    OldFileFmt.WriteFooter(wr, PermDAGFileType, PermDAGFileVersion);
  END Dump;

PROCEDURE DumpCompr(permdag: T; wr: Wr.T) RAISES {Huffman.OverFlow} =
  VAR
    code         : Code.T := Code.NewCode();
    firststate   : State := 1;
    laststate    : State;
    tail         : TEXT := "";
    f            : Huffman.FrequencyTable;
    h            : Huffman.T;
  BEGIN
    Wr.PutText(wr, ComprHeader); FPut.EOL(wr);
    code.InitWrite(wr);
    code.Write(WordSizeSize, permdag.format);
    WITH  com = permdag.comment,
          len = Text.Length(com)
          DO
      code.Write(WordSize, len);
      FOR i := 0 TO len - 1 DO
        code.Write(CharSize, ORD(Text.GetChar(com, i)))
      END
    END;
    WITH dag      = permdag.dag,
         m        = dag.MaxState(),
         msize    = Code.BitSize(m),
         permlast = LAST(permdag.charperm^),
         rdsize   = Code.BitSize(permlast),
         nulrd    = permdag.chartolettertable[FinalBitChar]
        DO
      code.Write(CharSize, rdsize);
      code.Write(rdsize, permlast);
      FOR c := FIRST(CHAR) TO LAST(CHAR) DO
        WITH t = permdag.chartolettertable[c],
             codet = t # InexistentLetter
             DO
          code.Write(BooleanSize, ORD(codet));
          IF codet THEN
            code.Write(rdsize, t)
          END
        END
      END;
      code.Write(WordSizeSize, msize - 1);
      WITH root = permdag.root DO
        code.Write(WordSizeSize, NUMBER(root^));
        FOR i := 0 TO LAST(root^) DO
          code.Write(msize, root[i])
        END
      END;

      f := NEW(Huffman.FrequencyTable, 1 + permlast);
      FOR c := 0 TO permlast DO f[c] := 0 END;
      FOR s := 1 TO m DO INC(f[dag.Last(s).rd]) END;
      h := Huffman.Tbl2Huffman(f);
      h.Dump(code);
      
      code.Write(msize, m);
      laststate := MIN(m, 2);
      FOR statesize := 0 TO msize DO
        FOR s := firststate TO laststate DO
          WITH arc      = dag.Last(s),
               rd       = 0 + arc.rd,
               dest     = arc.dest,
               class    = ClassifyAddress(s, dest),
               codedest = class = AddressClass.Jmp
               DO
            <* ASSERT (rd = nulrd) = (class = AddressClass.Nul) *>
            h.Compress(VAL(rd, CHAR), code); (* code.Write(rdsize, rd); *)
            code.Write(BooleanSize, ORD(codedest));
            IF codedest THEN
              code.Write(statesize, dest - 1)
            END
          END;
          WITH rest     = dag.Rest(s),
               class    = ClassifyAddress(s, rest),
               seq      = class = AddressClass.Seq,
               coderest = class = AddressClass.Jmp
               DO
            code.Write(BooleanSize, ORD(seq));
            IF NOT seq THEN
              code.Write(BooleanSize, ORD(coderest));
              IF coderest THEN
                code.Write(statesize, rest - 1)
              END
            END
          END
        END;
        IF laststate = m THEN
          code.Fin(tail);
          EXIT
        ELSE
          firststate := 1 + laststate;
          laststate := MIN(m, 2 * laststate - 1)
        END
      END
    END;
    Wr.PutText(wr, ComprTrailer); FPut.EOL(wr);
  END DumpCompr;
  
PROCEDURE DumpFixedBin(permdag: T; wr: Wr.T) =
  VAR
    code         : Code.T := Code.NewCode();
  BEGIN
    Wr.PutText(wr, ComprHeader); FPut.EOL(wr);
    code.InitWrite(wr);
    code.Write(WordSizeSize, permdag.format);
    WITH  com = permdag.comment,
          len = Text.Length(com)
          DO
      code.Write(WordSize, len);
      FOR i := 0 TO len - 1 DO
        code.Write(CharSize, ORD(Text.GetChar(com, i)))
      END
    END;
    WITH dag      = permdag.dag,
         m        = dag.MaxState(),
         msize    = Code.BitSize(m),
         permlast = LAST(permdag.charperm^),
         rdsize   = Code.BitSize(permlast)
         (* , nulrd    = permdag.chartolettertable[FinalBitChar] *)
        DO
      code.Write(CharSize, rdsize);
      code.Write(rdsize, permlast);
      FOR c := FIRST(CHAR) TO LAST(CHAR) DO
        WITH t = permdag.chartolettertable[c],
             codet = t # InexistentLetter
             DO
          code.Write(BooleanSize, ORD(codet));
          IF codet THEN
            code.Write(rdsize, t)
          END
        END
      END;
      code.Write(WordSizeSize, msize - 1);
      WITH root = permdag.root DO
        code.Write(WordSizeSize, NUMBER(root^));
        FOR i := 0 TO LAST(root^) DO
          code.Write(msize, root[i])
        END
      END;

      code.Write(msize, m);
      FOR s := 1 TO m DO
        WITH arc      = dag.Last(s)
             DO
          code.Write(rdsize, arc.rd);
          code.Write(msize,  arc.dest);
          code.Write(msize,  dag.Rest(s))
        END
      END
    END;
    Wr.PutText(wr, ComprTrailer); FPut.EOL(wr);
  END DumpFixedBin;
  
PROCEDURE ExpandState(    permdag : T; 
                            s       : State;
                        VAR final   : BOOLEAN;
                        VAR chars   : REF ARRAY OF CHAR) =
  VAR
    fullchars : ARRAY CHAR OF CHAR;
    nchars    : CHAR := FIRST(CHAR);
  BEGIN
    final := FALSE;
    WITH dag  = permdag.dag,
         perm = permdag.charperm
         DO
      WHILE s > 0 DO
        WITH arc = dag.Last(s),
             c   = perm[arc.rd]
             DO
          IF c = FinalBitChar AND arc.dest = 0 THEN
            final := TRUE
          ELSE
            fullchars[nchars] := c;
            INC(nchars)
          END
        END;
        s := dag.Rest(s)
      END
    END;
    WITH n = ORD(nchars) DO
      chars  := NEW(REF ARRAY OF CHAR, n);
      chars^ := SUBARRAY(fullchars, 0, n)
    END
  END ExpandState;
  
PROCEDURE Fold(permdag     : T;
                 backupname  : TEXT;
                 backuplevel : CARDINAL;
                 restore     : BOOLEAN
                ) RAISES {Full, Huffman.OverFlow} =
  TYPE
    ExtendedState    = [-1..LAST(State)];
    StateTable       = ARRAY OF ExtendedState;
    StateDescr       = RECORD
        outdegree    : CARDINAL := 0;
        descr        : REF StateTable;
      END;
    CardTable        = ARRAY OF CARDINAL;
  CONST
    flag             : ExtendedState  = -1;
  VAR
    alarmtime    : Time.T;
    bottomdescr  : StateDescr;
    count        : REF StateTable;
    currenttime  : Time.T;
    dagnumber    : CARDINAL;
    dagsize      : CARDINAL;
    equiv        : REF ARRAY OF CARDINAL;
    extradag     : DAG.T;
    indegree     : REF ARRAY OF CARDINAL;
    initdescr    : REF StateTable;
    initialsize  : CARDINAL;
    lasts        : State;
    matches      : REF StateTable;
    permlast     : CARDINAL;
    permlength   : CARDINAL;
    rdstack      : REF CardTable;
    statedescr   : StateDescr;
    workpermdag  : T := New();
    zeroroot     : REF ARRAY OF State := NEW(REF ARRAY OF State, 0);
  
  PROCEDURE Visit(READONLY stovisit: State) =
    VAR
      s              : State;
      imax           : CARDINAL;
    BEGIN
      statedescr.outdegree  := 0; statedescr.descr^  := initdescr^;
      bottomdescr.outdegree := 0; bottomdescr.descr^ := initdescr^;
      WITH e       = equiv[stovisit], 
           dag     = permdag.dag,
           workdag = workpermdag.dag,
           staout  = statedescr.outdegree,
           botout  = bottomdescr.outdegree
           DO
        s := stovisit;
        REPEAT
          INC(staout);
          WITH arc = dag.Last(s) DO
            statedescr.descr[arc.rd] := equiv[arc.dest];
            rdstack[permlength - staout] := arc.rd
          END;
          s := dag.Rest(s)
        UNTIL s = 0;
        SUBARRAY(rdstack^, 0, staout) := 
          SUBARRAY(rdstack^, permlength - staout, staout);
        FOR t := 1 TO workdag.MaxState() DO
          WITH mr  = matches[workdag.Rest(t)],
               m   = matches[t],
               arc = workdag.Last(t)
               DO
            IF mr = flag OR statedescr.descr[arc.rd] # arc.dest THEN
              m := flag
            ELSE
              m := 1 + mr;
              IF botout < m THEN
                botout := m;
                e := t;
                IF botout = staout THEN RETURN END
              END
            END
          END
        END;
        IF botout = 0 THEN e := 0 END;
        s := e;
        WHILE s # 0 DO
          WITH rd = workdag.Last(s).rd + 0, d = statedescr.descr[rd] DO
            bottomdescr.descr[rd] := d;
            d := flag;
            DEC(staout);
            FOR i := 0 TO staout DO
              IF rdstack[i] = rd THEN
                SUBARRAY(rdstack^, i, staout - i) := 
                  SUBARRAY(rdstack^, 1 + i, staout - i);
                EXIT
              END
            END
          END;
          s := workdag.Rest(s)
        END;
        WHILE staout > 1 DO
          count^ := initdescr^;
          FOR t := 1 TO dag.MaxState() DO
            WITH m   = matches[t],
                 arc = dag.Last(t)
                 DO
              m := matches[dag.Rest(t)];
              IF bottomdescr.descr[arc.rd] = equiv[arc.dest] THEN
                INC(m)
              END;
              IF m = botout AND indegree[t] > 1 AND t > stovisit THEN
                s := t;
                REPEAT
                  WITH arc = dag.Last(s) DO
                    IF statedescr.descr[arc.rd] = equiv[arc.dest] THEN
                      INC(count[arc.rd])
                    END
                  END;
                  s := dag.Rest(s)
                UNTIL s = 0
              END
            END
          END;
          imax := 0;
          FOR i := 1 TO staout DO
            IF count[rdstack[i]] > count[rdstack[imax]] THEN imax := i END
          END;
          IF count[rdstack[imax]] = flag THEN EXIT END;
          WITH cmax = rdstack[imax],
               d    = statedescr.descr[cmax]
               DO
            TRY
              e := workdag.Append(Arc{rd := cmax, dest := d}, e)
            EXCEPT
              Full => <* ASSERT FALSE *>
            END;
            bottomdescr.descr[cmax] := d;
            d := flag
          END;
          INC(botout);
          DEC(staout);
          SUBARRAY(rdstack^, imax, staout - imax) :=
            SUBARRAY(rdstack^, 1 + imax, staout - imax)
        END;
        FOR i := 0 TO staout - 1 DO
          WITH r = rdstack[i] DO
            TRY
              e := workdag.Append(Arc{rd := r, dest := statedescr.descr[r]}, e)
            EXCEPT
              Full => <* ASSERT FALSE *>
            END
          END
        END
      END
    END Visit;

  BEGIN
    currenttime := Time.Now();
    IF restore THEN
      FoldRestore(permdag, initialsize, lasts, equiv, workpermdag, backupname,
              backuplevel);
      alarmtime    := currenttime + BackupPeriod
    ELSE
      initialsize  := permdag.dag.MaxState();
      lasts        := 0;
      equiv        := NEW(REF ARRAY OF CARDINAL, 1 + initialsize);
      equiv[0]     := 0; 
      workpermdag.textperm           := permdag.textperm;
      workpermdag.charperm           := permdag.charperm;
      workpermdag.chartolettertable  := permdag.chartolettertable;
      alarmtime    := currenttime
    END;
    dagsize      := permdag.dag.MaxState();
    dagnumber    := 1 + dagsize;
    workpermdag.dag.Expand(dagnumber);
    indegree     := NEW(REF ARRAY OF CARDINAL, dagnumber);
    matches      := NEW(REF StateTable, dagnumber);
    matches[0]   := 0;
    
    permlast     := LAST(permdag.charperm^);
    permlength   := 1 + permlast;
    initdescr    := NEW(REF StateTable, permlength);
    FOR i := 0 TO permlast DO initdescr[i] := flag END; 
    bottomdescr  := StateDescr{
      outdegree := 0, descr := NEW(REF StateTable, permlength)};
    statedescr   := StateDescr{
      outdegree := 0, descr := NEW(REF StateTable, permlength)};
    count        := NEW(REF StateTable, permlength);
    rdstack      := NEW(REF CardTable,  permlength);
    REPEAT
      dagsize := permdag.dag.MaxState();
      Wr.PutText(stderr, Fmt.Int(dagsize) & "\n");
      Wr.Flush(stderr);
      indegree[0] := 0;
      FOR s := 1 TO dagsize DO
        indegree[s] := 0;
        INC(indegree[permdag.dag.Last(s).dest], 2);
        INC(indegree[permdag.dag.Rest(s)])
      END;
      FOR i := 0 TO LAST(permdag.root^) DO
        INC(indegree[permdag.root[i]], 2)
      END;
      FOR s := 1 + lasts TO dagsize DO
        equiv[s] := 0;
        IF indegree[s] > 1 THEN
          Visit(s);
          currenttime := Time.Now();
          IF currenttime >= alarmtime THEN
            FoldBackUp(permdag, initialsize, s, equiv, workpermdag,
                       backupname, backuplevel);
            alarmtime := currenttime + BackupPeriod
          END
        END
      END;
      WITH root = permdag.root DO
        FOR i := 0 TO LAST(root^) DO
          WITH r = root[i] DO
            r := equiv[r]
          END
        END
      END;
      permdag.dag.Crunch(zeroroot^);
      extradag := permdag.dag; 
      permdag.dag := workpermdag.dag;
      workpermdag.dag := extradag;
      lasts := 0
    UNTIL dagsize = permdag.dag.MaxState();
    Wr.PutText(stderr, "\ninitial = " & Fmt.Int(initialsize) &
      " final = " & Fmt.Int(dagsize) & " => " & 
      Fmt.Int(ROUND(FLOAT(dagsize) * 100.0 / FLOAT(initialsize))) &
      "%\n");
    Wr.Flush(stderr)
  END Fold;
  
PROCEDURE FoldBackUp(    permdag     : T;
                           initialsize : CARDINAL;
                           lasts       : State;
                           equiv       : REF ARRAY OF CARDINAL;
                           workpermdag : T;
                           backupname  : TEXT;
                       VAR backuplevel : CARDINAL
                      ) RAISES {Huffman.OverFlow}=
  VAR
    date      : Date.T;
    code      : Code.T   := Code.NewCode();
    tail      : TEXT;
  BEGIN
    date := Util.GetDate();
    Wr.PutText(stderr, Util.FmtDate(date) & " backup = " &
      Fmt.Int(backuplevel) & ": " & 
      Fmt.Int(ROUND(100.0 * FLOAT(lasts) / FLOAT(permdag.dag.MaxState()))) & 
      "%");
    Wr.Flush(stderr);
    WITH 
      size = Code.BitSize(permdag.dag.MaxState()),
      backupwriter = FileWr.Open( 
      backupname & Fmt.Int(backuplevel)) 
    DO
      permdag.DumpCompr(backupwriter);
      code.InitWrite(backupwriter);
      code.Write(WordSize, initialsize);
      code.Write(size, lasts);
      FOR s := 0 TO lasts DO
        code.Write(size, equiv[s])
      END;
      tail := "";
      code.Fin(tail);
      workpermdag.DumpCompr(backupwriter);
      Wr.Close(backupwriter);
      Wr.PutText(stderr, ": done\n");
      Wr.Flush(stderr)
    END;
    backuplevel := 1 - backuplevel
  END FoldBackUp;
  
PROCEDURE FoldRestore(    permdag     : T;
                        VAR initialsize : CARDINAL;
                        VAR lasts       : State;
                        VAR equiv       : REF ARRAY OF CARDINAL;
                            workpermdag : T;
                            backupname  : TEXT;
                        VAR backuplevel : CARDINAL
                       ) RAISES {Full} =
  VAR
    code : Code.T   := Code.NewCode();
    tail : TEXT;
  <* FATAL Rd.EndOfFile *>
  BEGIN
    WITH backupreader = FileRd.Open(
           backupname & Fmt.Int(backuplevel)) DO
      Copy(InCodeCompr(backupreader), permdag);
      code.InitRead(backupreader);
      initialsize := code.Read(WordSize);
      WITH size       = permdag.dag.MaxState(),
           sizelength = Code.BitSize(size) DO
        lasts := code.Read(sizelength);
        equiv := NEW(REF ARRAY OF CARDINAL, 1 + size);
        FOR s := 0 TO lasts DO
          equiv[s] := code.Read(sizelength)
        END
      END;
      code.Fin(tail);
      Copy(InCodeCompr(backupreader), workpermdag);
      Rd.Close(backupreader)
    END;
    backuplevel :=  1 - backuplevel
  END FoldRestore;
  
PROCEDURE InCodeCompr(rd: Rd.T): T 
    RAISES {Full} =
  VAR
    permdag      : T := New();
    firststate   : State := 1;
    laststate    : State;
    dest         : State;
    rest         : State;
    
    code         : Code.T := Code.NewCode();
    rdsize       : Letter;
    msize        : CARDINAL;
    tail         : TEXT;
    h            : Huffman.T;
    rdVar        : CARDINAL;
  <* FATAL Rd.EndOfFile *>
  BEGIN
    FGet.Match(rd, ComprHeader); FGet.EOL(rd);
    code.InitRead(rd);
    WITH format = permdag.format,
         com    = permdag.comment DO
      format := code.Read(WordSizeSize);
      IF format >= 1 THEN
        com := "";
        FOR i := 1 TO code.Read(WordSize) DO
          com := com & Text.FromChar(VAL(code.Read(CharSize), CHAR))
        END
      ELSE
        com := DefaultCommentText
      END
    END;
    rdsize := code.Read(CharSize);
    WITH charperm   = permdag.charperm DO
      charperm := NEW(REF ARRAY OF CHAR, 1 + code.Read(rdsize));
      FOR c := FIRST(CHAR) TO LAST(CHAR) DO
        WITH codet = VAL(code.Read(BooleanSize), BOOLEAN),
             t     = permdag.chartolettertable[c]
             DO
          IF codet THEN
            t := code.Read(rdsize);
            charperm[t] := c
          ELSE
            t := InexistentLetter
          END
        END
      END;
      permdag.textperm := Text.FromChars(charperm^)
    END;   
    msize := 1 + code.Read(WordSizeSize);
    WITH root = permdag.root DO
      root := NEW(REF ARRAY OF State, code.Read(WordSizeSize));
      FOR i := 0 TO LAST(root^) DO
        root[i] := code.Read(msize)
      END
    END;

    IF permdag.format >= 4 THEN 
      h := Huffman.Load(code, NUMBER(permdag.charperm^));
    END;

    WITH m     = code.Read(msize),
         nulrd = permdag.chartolettertable[FinalBitChar],
         dag   = permdag.dag 
         DO
      dag := DAG.New(m);
      laststate := MIN(m, 2);
      FOR statesize := 0 TO msize DO
        FOR s := firststate TO laststate DO
          IF permdag.format >= 4 
          THEN rdVar := ORD(h.Uncompress(code))
          ELSE rdVar := code.Read(rdsize)
          END;
          IF VAL(code.Read(BooleanSize), BOOLEAN) THEN
            dest := 1 + code.Read(statesize)
          ELSIF rdVar = nulrd THEN
            dest := 0
          ELSE
            dest := s - 1
          END;
          IF VAL(code.Read(BooleanSize), BOOLEAN) THEN
            rest := s - 1
          ELSIF VAL(code.Read(BooleanSize), BOOLEAN) THEN
            rest := 1 + code.Read(statesize)
          ELSE
            rest := 0
          END;
          EVAL dag.Append(Arc{rd := rdVar, dest := dest}, rest)
        END;
        IF laststate = m THEN
          code.Fin(tail);
          EXIT
        ELSE
          firststate := 1 + laststate;
          laststate := MIN(m, 2 * laststate - 1)
        END
      END
    END;
    IF permdag.format >= 2 THEN
      IF permdag.format >= 3 THEN
        tail := ""
      END;
      
    END;
    permdag.format := PermDAGFormat;

    (* Check trailer: *)
    TRY
      WITH line = Rd.GetLine(rd) DO tail := tail & line END
    EXCEPT 
      Rd.EndOfFile => (*OK*)
    END;
    <* ASSERT Text.Equal(tail, ComprTrailer) *>

    RETURN permdag
  END InCodeCompr;
  
PROCEDURE IncRoots(permdag: T) =
  BEGIN
    WITH root    = permdag.root,
         n       = NUMBER(root^),
         newroot = NEW(REF ARRAY OF State, 1 + n)
         DO
      SUBARRAY(newroot^, 0, n) := root^;
      root := newroot;
      root[n] := 0
    END
  END IncRoots;
  
PROCEDURE IsFinal(permdag: T; s: State): BOOLEAN =
  VAR
    final : BOOLEAN;
    chars : REF ARRAY OF CHAR;
  BEGIN
    ExpandState(permdag, s, final, chars);
    RETURN final
  END IsFinal;
  
PROCEDURE Load(rd: Rd.T): T =
  <* FATAL Rd.EndOfFile *>
  VAR permdag: T := New();
  BEGIN
    OldFileFmt.ReadHeader(rd, PermDAGFileType, PermDAGFileVersion);
    
    permdag.format := NGet.Int(rd, FormatPrefix);
    FGet.EOL(rd);
    
    permdag.comment := OldFileFmt.ReadComment(rd, '|');
    
    (* Careful - blanks are significant here! *)
    FGet.Match(rd, PermPrefix); 
    FGet.Match(rd, " = ");
    permdag.textperm := Rd.GetLine(rd);
    BuildPermDagTables(permdag);
    
    WITH 
      nRoots = NGet.Int(rd, RootPrefix) + 1,
      r = NEW(REF ARRAY OF State, nRoots)
    DO
      FGet.Colon(rd);
      FOR i := 0 TO nRoots-1 DO r[i] := FGet.Int(rd) END;
      permdag.root := r
    END;
    FGet.EOL(rd);
    
    permdag.dag := DAG.Load(rd, 0);
    
    OldFileFmt.ReadFooter(rd, PermDAGFileType, PermDAGFileVersion);
    permdag.format := PermDAGFormat;
    RETURN permdag
  END Load;
  
PROCEDURE LoadCompr(rd: Rd.T): T RAISES {Full} =
  BEGIN
    RETURN InCodeCompr(rd)
  END LoadCompr;
  
PROCEDURE MakeTransition(permdag: T; s: State; c: CHAR): State =
  BEGIN
    WITH dag    = permdag.dag,
         letter = permdag.chartolettertable[c]
         DO
      WHILE s > 0 DO
        WITH arc = dag.Last(s) DO
          IF arc.rd = letter THEN
            RETURN arc.dest
          END
        END;
        s := dag.Rest(s)
      END
    END;
    RETURN s
  END MakeTransition;
  
PROCEDURE New(): T =
  BEGIN
    WITH permdag = NEW(T,
      format            := PermDAGFormat,
      comment           := DefaultCommentText,
      textperm          := Text.FromChar(FinalBitChar),
      charperm          := NEW(REF ARRAY OF CHAR, 1),
      chartolettertable := ARRAY CHAR OF ExtendedLetter{InexistentLetter, ..},
      root              := NEW(REF ARRAY OF State, 1),
      dag               := DAG.New(0))
    DO
      permdag.charperm[0] := FinalBitChar;
      permdag.chartolettertable [FinalBitChar] := 0;
      permdag.root[0] := 0;
      RETURN permdag
    END
  END New;
  
EXCEPTION BadPermDAG;
  
PROCEDURE Prm2Red(permdag: T; wr: Wr.T) RAISES {Full} =
  VAR letter  : CARDINAL;
  <* FATAL BadPermDAG *>
  BEGIN
    OldFileFmt.WriteHeader(wr, ReducedFileType, ReducedFileVersion);
    
    OldFileFmt.WriteComment(wr, permdag.comment, '|');
    
    WITH dag    = permdag.dag,
         perm   = permdag.charperm,
         root   = permdag.root,
         m      = permdag.dag.MaxState(),
         reddag = DAG.New(m)
         DO
      <* ASSERT LAST(root^) = 0 *>
      FOR s := 1 TO m DO
        WITH arc  = dag.Last(s),
             c    = perm[arc.rd],
             dest = arc.dest,
             rest = dag.Rest(s)
             DO
          IF c = FinalBitChar AND dest = 0 THEN
            IF rest = 0 THEN
              letter := 0
            ELSE
              RAISE BadPermDAG
            END
          ELSE
            letter := ORD(c)
          END;
          EVAL reddag.Append(Arc{rd := letter, dest := dest}, rest);
          IF rest # 0 AND reddag.Last(rest).rd >= letter THEN
            RAISE BadPermDAG
          END;
        END
      END;
      
      NPut.Int(wr, "root", root[0]); FPut.EOL(wr);
      DAG.Dump (wr, reddag);
    END;
    OldFileFmt.WriteFooter(wr, ReducedFileType, ReducedFileVersion);
  END Prm2Red;
  
PROCEDURE Red2Prm(rd: Rd.T): T RAISES {Full} =
  VAR
    permdag   : T := New();
    usedchars : SET OF CHAR := SET OF CHAR{};
    equivchar : CHAR;
  BEGIN
    OldFileFmt.ReadHeader(rd, ReducedFileType, ReducedFileVersion);
    
    permdag.comment := OldFileFmt.ReadComment(rd, '|');
    
    WITH 
      redRoot = NGet.Int(rd, "root"),
      r = NEW(REF ARRAY OF State, 1)
    DO
      r[0] := redRoot;
      permdag.root := r
    END;
    FGet.EOL(rd);
    
    WITH dag    = permdag.dag,
         t      = permdag.textperm,
         reddag = DAG.Load(rd),
         m      = reddag.MaxState()
         DO
      dag := DAG.New(m);
      FOR s := 1 TO m DO
        WITH rd = reddag.Last(s).rd + 0 DO
          IF rd # 0 THEN
            usedchars := usedchars + SET OF CHAR{VAL(rd, CHAR)}
          END
        END
      END;
      <* ASSERT NOT (FinalBitChar IN usedchars) *>
      FOR c := FIRST(CHAR) TO LAST(CHAR) DO
        IF c IN usedchars THEN
          t := t & Text.FromChar(c)
        END
      END;
      BuildPermDagTables(permdag);
      FOR s := 1 TO m DO
        WITH arc = reddag.Last(s),
             rd  = arc.rd + 0
             DO
          IF rd = 0 THEN
            equivchar := FinalBitChar
          ELSE
            equivchar := VAL(arc.rd, CHAR)
          END;
          EVAL dag.Append(
            Arc{rd   := permdag.chartolettertable[equivchar],
                dest := arc.dest}, 
            reddag.Rest(s))
        END
      END
    END;
    OldFileFmt.ReadFooter(rd, ReducedFileType, ReducedFileVersion);
    RETURN permdag
  END Red2Prm;

PROCEDURE Root(permdag: T; level : CARDINAL := 0): State =
  BEGIN
    WITH root = permdag.root DO
      <* ASSERT level <= LAST(root^) *>
      RETURN root[level]
    END
  END Root;
  
PROCEDURE SortPrm(permdag: T) RAISES {Full} =
  VAR
    Present       : ARRAY CHAR OF BOOLEAN := ARRAY CHAR OF BOOLEAN {FALSE, ..};
    WorkPermDag   : T := New();
  BEGIN
    WITH dag      = permdag.dag,
         m        = dag.MaxState(),
         charperm = permdag.charperm,
         workdag  = WorkPermDag.dag
    DO
      WorkPermDag.comment := permdag.comment;
      WorkPermDag.root    := NEW(REF ARRAY OF State, NUMBER(permdag.root^));
      WorkPermDag.root^   := permdag.root^;
      workdag.Expand(1 + m);
      FOR s := 1 TO m DO
        Present[charperm[dag.Last(s).rd]] := TRUE
      END;
      Present[FinalBitChar] := FALSE;
      WITH t = WorkPermDag.textperm DO
        FOR c := FIRST(CHAR) TO LAST(CHAR) DO
          IF  Present[c] THEN t := t & Text.FromChar(c) END
        END
      END;
      BuildPermDagTables(WorkPermDag);
      WITH table = WorkPermDag.chartolettertable DO
        TRY
          FOR s := 1 TO m DO
            WITH arc = dag.Last(s), r = arc.rd + 0, d = arc.dest DO
              EVAL workdag.Append(Arc{rd := table[charperm[r]], dest := d}, 
                                  dag.Rest(s))
            END
          END
        EXCEPT
          Full => <* ASSERT FALSE *>
        END
      END
    END;
    Copy(WorkPermDag, permdag)
  END SortPrm;
  
PROCEDURE Spell(permdag: T; wr: Wr.T; level: CARDINAL := 0) =
  VAR
    Flag          : State := 1 + permdag.dag.MaxState();
    Word          : ARRAY [1..100] OF CHAR;
    WordLength    : CARDINAL := 0;

  PROCEDURE Visit(s: State) =
    VAR
      StateDescr : ARRAY CHAR OF State := ARRAY CHAR OF State{Flag, ..};
    BEGIN
      WITH dag = permdag.dag,
           perm = permdag.charperm
           DO
        WHILE s # 0 DO
          WITH arc = dag.Last(s), c = perm[arc.rd] DO
            IF c = FinalBitChar THEN
              Wr.PutString(wr, SUBARRAY(Word, 0, WordLength));
              Wr.PutChar(wr, '\n')
            ELSE
              StateDescr[c] := arc.dest
            END
          END;
          s := dag.Rest(s)
        END;
        INC(WordLength);
        WITH w = Word[WordLength] DO
          FOR c := FIRST(CHAR) TO LAST(CHAR) DO
            WITH dest = StateDescr[c] DO
              IF dest # Flag THEN
                w := c;
                Visit(dest)
              END
            END
          END
        END;
        DEC(WordLength)
      END
    END Visit;

  BEGIN
    WITH root = permdag.root DO
      <* ASSERT level <= LAST(root^) *>
      Visit(root[level])
    END
  END Spell;
  
PROCEDURE Statistics(permdag: T; wr: Wr.T) =
  TYPE
    TypeCounter = RECORD
        Visited  : BOOLEAN;
        IsState  : BOOLEAN;
        IsFinal  : BOOLEAN;
        Arcs     : CARDINAL;
        Letters  : CARDINAL;
        Words    : CARDINAL
      END;
  VAR
    States       : CARDINAL := 0;
    Finals       : CARDINAL := 0;
    Arcs         : CARDINAL := 0;
    Letters      : CARDINAL := 0;
    Words        : CARDINAL := 0;
    Counter      : REF ARRAY OF TypeCounter;
    dag          : DAG.T := permdag.dag;
    perm         : REF ARRAY OF CHAR := permdag.charperm;

  PROCEDURE Visit(s: State) =
    BEGIN
      WITH c = Counter[s] DO
        IF NOT c.Visited THEN
          c.Visited := TRUE;
          WITH arc = dag.Last(s),
               rd  = perm[arc.rd],
               dest = arc.dest,
               rest = dag.Rest(s)
               DO
            Visit(dest);
            Visit(rest);
            Counter[dest].IsState := TRUE;
            WITH d = Counter[dest],
                 r = Counter[rest]
              DO
              c.Arcs     := r.Arcs;
              c.Letters  := d.Letters + r.Letters;
              c.Words    := d.Words + r.Words;
              IF rd = FinalBitChar AND dest = 0 THEN
                c.IsFinal := TRUE;
                INC(c.Words)
              ELSE
                INC(c.Arcs);
                c.IsFinal := Counter[rest].IsFinal;
                INC(c.Letters, d.Words)
              END
            END
          END
        END
      END
    END Visit;

  BEGIN
    WITH m    = dag.MaxState(),
         root = permdag.root
      DO
      Counter := NEW(REF ARRAY OF TypeCounter, 1 + m);
      Counter[0] := TypeCounter{Visited := TRUE, IsState := FALSE,
        IsFinal := FALSE, Arcs := 0, Letters := 0, Words   := 0};
      FOR i := 1 TO m DO
        Counter[i].Visited := FALSE
      END;
      FOR i := 0 TO LAST(root^) DO
        WITH r = root[i] DO
          Counter[r].IsState := TRUE;
          Visit(r)
        END
      END;
      FOR i := 1 TO m DO
        WITH c = Counter[i] DO
          IF c.IsState THEN
            INC(States);
            IF c.IsFinal THEN INC(Finals) END;
            INC(Arcs, c.Arcs)
          END
        END
      END;
      FOR i := 0 TO LAST(root^) DO
        WITH c = Counter[root[i]] DO
          INC(Letters, c.Letters);
          INC(Words, c.Words)
        END
      END;
      Wr.PutText(wr,
        Fmt.F("%9s %9s %9s %9s", "strings", "letters", "states", "finals") &
        Fmt.F("%9s %9s %9s", "arcs", "sub-sts", "lets/arc") &"\n" &
        Fmt.F("%9s %9s %9s %9s",
          "--------", "--------", "--------", "--------") &
        Fmt.F("%9s %9s %9s", "--------", "--------", "--------") &"\n" &
        Fmt.F("%9s %9s %9s %9s",
          Fmt.Int(Words), Fmt.Int(Letters), Fmt.Int(States), Fmt.Int(Finals)) &
        Fmt.F("%9s %9s %9s", Fmt.Int(Arcs), Fmt.Int(m),
          Fmt.Real(FLOAT(Letters) / FLOAT(Arcs),
            prec := 3, style := Fmt.Style.Fix))
        & "\n");
      Wr.Flush(wr)
    END
  END Statistics;
  
PROCEDURE UnFold(permdag: T) RAISES {Full} =
  VAR
    dagsize : State              := permdag.dag.MaxState();
    newdag  : DAG.T              := DAG.New(dagsize);
    equiv   : REF ARRAY OF State := NEW(REF ARRAY OF State, 1 + dagsize);
    
    PROCEDURE Visit(s: State): State RAISES {Full} =
    TYPE
      CharIndices = [ORD(FIRST(CHAR))..ORD(LAST(CHAR))];
    VAR
      UsedChars : SET OF CharIndices          := SET OF CharIndices{};
      Dests     : ARRAY CharIndices OF State;
    BEGIN
      WITH e    = equiv[s],
           dag  = permdag.dag
           DO
        IF e # 0 OR s = 0 THEN RETURN e END;
        REPEAT
          WITH arc = dag.Last(s) DO
            UsedChars := UsedChars + SET OF CharIndices{arc.rd};
            Dests[arc.rd]  := arc.dest
          END;
          s := dag.Rest(s)
        UNTIL s = 0;
        FOR c := FIRST(CharIndices) TO LAST(CharIndices) DO
          IF c IN UsedChars THEN
            LOOP
              TRY
                e := newdag.Append(Arc{rd := c, dest := Visit(Dests[c])}, e);
                EXIT
              EXCEPT
                Full => newdag.Expand(2 + 11 * newdag.MaxAllocState() DIV 10)
              END
            END
          END
        END;
        RETURN e
      END
    END Visit;
  
  BEGIN
    FOR s := 0 TO dagsize DO
      equiv[s] := 0
    END;
    FOR i := 0 TO LAST(permdag.root^) DO
      WITH r = permdag.root[i] DO
        r := Visit(r)
      END
    END;
    permdag.dag := newdag
  END UnFold;

BEGIN
END PermDAG.
