MODULE Huffman;

IMPORT Rd, Code, Word;

TYPE
  Bit = [0..1];
  Tree = REF TreeNode;
  TreeNode       = RECORD
      s: ARRAY Bit OF Tree := ARRAY Bit OF Tree{NIL, NIL};
      c: CHAR;
    END;
  CompressionTableEntry = RECORD
      l: Code.TypeSize;
      n: Word.T;
    END;
  CodedTreeRecord       = RECORD
      a: REF ARRAY OF BOOLEAN;
      l: REF ARRAY OF CHAR;
    END;
REVEAL
  T = Public BRANDED OBJECT
      TableSize        : CARDINAL;
      CharSize         : CARDINAL;
      Root             : Tree;
      CompressionTable : REF ARRAY OF CompressionTableEntry;
      TreeSize         : CARDINAL;
      CodedTree        : CodedTreeRecord
    OVERRIDES
      Compress         := Compress;
      Uncompress       := Uncompress;
      Dump             := Dump
    END;
  
VAR
  BooleanSize := Code.BooleanSize();
  
  PROCEDURE New(TblSiz: CARDINAL): T =
  BEGIN
    WITH h = NEW(T) DO
      h.TableSize        := TblSiz;
      h.CharSize         := Code.BitSize(h.TableSize - 1);
      h.TreeSize         := 1 + 2 * h.TableSize;
      h.Root             := NIL;
      h.CompressionTable := NEW(REF ARRAY OF CompressionTableEntry, h.TableSize);
      h.CodedTree        := CodedTreeRecord
        {a := NEW(REF ARRAY OF BOOLEAN, h.TreeSize),
         l := NEW(REF ARRAY OF CHAR,  h.TableSize)};
      RETURN h
    END
  END New;

  PROCEDURE BuildTable(h: T) RAISES {OverFlow} =

    PROCEDURE AuxBuild(t: Tree; ct: CompressionTableEntry) RAISES 
      {OverFlow} =
      BEGIN
        WITH e = t.s[0] DO
          IF e = NIL THEN
            h.CompressionTable[ORD(t.c)] := ct
          ELSIF ct.l = Word.Size THEN
            RAISE OverFlow
          ELSE
            INC(ct.l); ct.n := Word.Shift(ct.n, 1); AuxBuild(e, ct);
            ct.n := Word.Or(ct.n, 1); AuxBuild(t.s[1], ct)
          END
        END;
      END AuxBuild;

  BEGIN
    AuxBuild(h.Root, CompressionTableEntry{0, 0});
  END BuildTable;

  PROCEDURE EncodeTree(h: T) RAISES {} =
  VAR
    na, nl : CARDINAL := 0;
    
    PROCEDURE AuxEncode(t: Tree) =
    BEGIN
      WITH e = t.s[0],
           x = e = NIL 
           DO
        h.CodedTree.a[na]  := x; INC(na);
        IF x THEN
          h.CodedTree.l[nl] := t.c; INC(nl)
        ELSE
          AuxEncode(e); AuxEncode(t.s[1])
        END
      END
    END AuxEncode;
    
  BEGIN
    AuxEncode(h.Root)
  END EncodeTree;

  PROCEDURE Tbl2Huffman (f: FrequencyTable): T RAISES {OverFlow} =
  VAR
    h: T := New(NUMBER(f^));
    FreqTab: REF ARRAY OF CARDINAL;
    NodeTab: REF ARRAY OF Tree;

    PROCEDURE Insert(n: CARDINAL) RAISES {} =
      VAR j: INTEGER;
      BEGIN
        <* ASSERT n > 0 *>
        VAR 
          fr := FreqTab[n];
          no := NodeTab[n];
        BEGIN
          j := 0;
          FOR i := n TO 1 BY -1 DO
            WITH fi = FreqTab[i-1] DO
              IF fi <= fr THEN 
                FreqTab[i] := fi;
                NodeTab[i] := NodeTab[i-1]
              ELSE 
                j := i; 
                EXIT
              END
            END
          END;
          FreqTab[j] := fr;
          NodeTab[j] := no;
        END;
      END Insert;

  BEGIN
    FreqTab := NEW(REF ARRAY OF CARDINAL, h.TableSize);
    NodeTab := NEW(REF ARRAY OF Tree,     h.TableSize);
    FOR n := 0 TO LAST(f^) DO
      FreqTab[n] := f[n];
      NodeTab[n] := NEW(Tree, c := VAL(n, CHAR));
    END;
    FOR n := 1 TO LAST(f^) DO Insert(n) END;
    FOR n := LAST(f^) TO 1 BY -1 DO
      Insert(n);
      INC(FreqTab[n-1], FreqTab[n]);
      WITH p = NodeTab[n-1],
           x = NEW(Tree) DO
        x^ := p^;
        p.s := ARRAY Bit OF Tree{x, NodeTab[n]}
      END
    END;
    h.Root := NodeTab[0];
    BuildTable(h);
    EncodeTree(h);
    RETURN h
  END Tbl2Huffman;
  
  PROCEDURE Load        (cd: Code.T; TblSiz: CARDINAL := NUMBER(CHAR)): T 
    RAISES {Rd.EndOfFile} =
  VAR
    h: T := New(TblSiz);
    na, nl : CARDINAL := 0;

    PROCEDURE AuxLoad(): Tree =
    VAR
      t: Tree := NEW(Tree);
      tr : BOOLEAN;
    BEGIN
      tr := h.CodedTree.a[na]; INC(na);
      IF tr THEN
        t.c := h.CodedTree.l[nl]; INC(nl)
      ELSE
        t.s[0] := AuxLoad(); t.s[1] := AuxLoad()
      END;
      RETURN t
    END AuxLoad;
    
  BEGIN
    FOR i := 0 TO h.TreeSize - 1 DO
      WITH a = h.CodedTree.a DO
        a[i] := VAL(cd.Read(BooleanSize), BOOLEAN)
      END
    END;
    FOR i := 0 TO h.TableSize - 1 DO
      WITH l = h.CodedTree.l DO
        l[i] := VAL(cd.Read(h.CharSize), CHAR);
      END
    END;
    h.Root := AuxLoad();
    RETURN h
  END Load;
  
  PROCEDURE Compress      (h: T; c: CHAR; cd: Code.T) =
  BEGIN
    WITH ct = h.CompressionTable[ORD(c)] DO
      cd.Write(ct.l, ct.n)
    END
  END Compress;
  
  PROCEDURE Uncompress    (h: T; cd: Code.T): CHAR RAISES {Rd.EndOfFile} =
  VAR
    t: Tree := h.Root;
  BEGIN
    WHILE t.s[0] # NIL  DO
      t := t.s[cd.Read(1)];
    END;
    RETURN t.c
  END Uncompress;
  
  PROCEDURE Dump          (h: T; cd: Code.T) =
  BEGIN
    FOR i := 0 TO h.TreeSize - 1 DO cd.Write(BooleanSize, ORD(h.CodedTree.a[i]))
    END; 
    FOR i := 0 TO h.TableSize - 1 DO cd.Write(h.CharSize, ORD(h.CodedTree.l[i]))
    END
  END Dump;

BEGIN
END Huffman.
