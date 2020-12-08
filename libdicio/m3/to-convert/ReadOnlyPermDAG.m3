MODULE ReadOnlyPermDAG;

IMPORT Code, Huffman;
IMPORT Rd, FGet, Thread, Text;

CONST
  CharSize     = BITSIZE(CHAR);
VAR
  BooleanSize  : CARDINAL := Code.BooleanSize();
  WordSize     : CARDINAL := Code.WordSize();
  WordSizeSize : CARDINAL := Code.WordSizeSize() - 1; (* excludes 0 *)
TYPE
  EntryRecord  = RECORD
      rd   : Letter;
      dest : State;
      rest : State
    END;
REVEAL
    T = Public BRANDED OBJECT
      comment            : TEXT;
      fields             : CARDINAL;
      textperm           : TEXT;
      charperm           : REF ARRAY OF CHAR;
      chartolettertable  : ARRAY CHAR OF ExtendedLetter;
      root               : REF ARRAY OF State;
      dag                : REF ARRAY OF EntryRecord;
    OVERRIDES
      Accepts            := Accepts;
      ExpandState        := ExpandState;
      IsFinal            := IsFinal;
      MakeTransition     := MakeTransition;
      Root               := Root
    END;

  PROCEDURE Accepts(permdag: T; s: State; word: TEXT): BOOLEAN =
  BEGIN
    WITH len = Text.Length(word),
         w   = word & Text.FromChar(FinalBitChar),
         dag = permdag.dag,
         t   = permdag.chartolettertable
         DO
      FOR i := 0 TO len DO
        WITH c = t[Text.GetChar(w, i)] DO
          LOOP
            IF s = 0 THEN RETURN FALSE END;
            WITH e = dag[s] DO
              IF e.rd = c THEN s := e.dest; EXIT
              ELSE s := e.rest
              END
            END
          END
        END
      END;
      RETURN s = 0
    END
  END Accepts;
  
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
        WITH e = dag[s],
             c   = perm[e.rd]
             DO
          IF c = FinalBitChar AND e.dest = 0 THEN
            final := TRUE
          ELSE
            fullchars[nchars] := c;
            INC(nchars)
          END;
          s := e.rest
        END
      END
    END;
    WITH n = ORD(nchars) DO
      chars  := NEW(REF ARRAY OF CHAR, n);
      chars^ := SUBARRAY(fullchars, 0, n)
    END
  END ExpandState;
  
  PROCEDURE IsFinal(permdag: T; s: State): BOOLEAN =
  VAR
    final : BOOLEAN;
    chars : REF ARRAY OF CHAR;
  BEGIN
    ExpandState(permdag, s, final, chars);
    RETURN final
  END IsFinal;
  
  PROCEDURE LoadCompr(rd: Rd.T): T =
  <* FATAL Rd.EndOfFile, Rd.Failure, Thread.Alerted *>
  VAR
    permdag      : T := NEW(T);
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
  BEGIN
    FGet.Match(rd, ComprHeader); FGet.EOL(rd);
    code.InitRead(rd);
    WITH fields = permdag.fields,
         com    = permdag.comment DO
      fields := code.Read(WordSizeSize);
      IF fields >= 1 THEN
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

    IF permdag.fields >= 4 THEN 
      h := Huffman.Load(code, NUMBER(permdag.charperm^));
    END;

    WITH m     = code.Read(msize),
         nulrd = permdag.chartolettertable[FinalBitChar],
         dag   = permdag.dag 
         DO
      dag := NEW(REF ARRAY OF EntryRecord, 1 + m);
      laststate := MIN(m, 2);
      FOR statesize := 0 TO msize DO
        FOR s := firststate TO laststate DO
          IF permdag.fields >= 4 
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
          dag[s] := EntryRecord{rd := rdVar, dest := dest, rest:= rest}
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
    IF permdag.fields >= 2 THEN
      IF permdag.fields >= 3 THEN
        tail := ""
      END;
    END;

    (* Check trailer: *)
    TRY
      WITH line = Rd.GetLine(rd) DO tail := tail & line END
    EXCEPT 
      Rd.EndOfFile => (*OK*)
    END;
    <* ASSERT Text.Equal(tail, ComprTrailer) *>

    RETURN permdag
  END LoadCompr;
  
  PROCEDURE MakeTransition(permdag: T; s: State; c: CHAR): State =
  BEGIN
    WITH dag    = permdag.dag,
         letter = permdag.chartolettertable[c]
         DO
      WHILE s > 0 DO
        WITH e = dag[s] DO
          IF e.rd = letter THEN RETURN e.dest END;
          s := e.rest
        END
      END
    END;
    RETURN 0
  END MakeTransition;
  
  PROCEDURE Root(permdag: T; level : CARDINAL := 0): State =
  BEGIN
    WITH root = permdag.root DO
      <* ASSERT level <= LAST(root^) *>
      RETURN root[level]
    END
  END Root;
  
BEGIN
END ReadOnlyPermDAG.
