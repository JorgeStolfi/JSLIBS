(* Oct-edge data structure. *)
(* Created 1993 by Jorge Stolfi and Rober M. Rosi  *)

MODULE Oct;

IMPORT Lex, Word, Wr, Rd, Fmt, FloatMode, Thread;
(* IMPORT Stdio; *)

TYPE 
  MarkBit = (* BITS 1 FOR *) BOOLEAN;

REVEAL 
  Edge = PublicEdge BRANDED
  OBJECT 
    marks: ARRAY [0..1] OF  MarkBit;  (* Used by EnumEdges, EnumVertices *)
    
    (* Element index [r] refers to arc a[r] = Arc{e,r}: *)
    enext: ARRAY [0..3] OF Edge;      (* Edge of Onext(a[r]) *)
    bnext: ARRAY [0..3] OF [0..7];    (* FRBits of Onext(a[r]) *)
  OVERRIDES  
    init := OctInit;
  END;
  (* 
    If "a = flip^f(rot^r(a0))", where "a0" is the reference arc
    of "a", "r" in [0..3], and "f" in [0..1], then "a.bits = (2r + f)".
    Therefore,

    | a.flip.bits = XOR(a.bits, 1)
    | If fbit is 0
    |   a.rot.bits = (a.bits + 2) MOD 8 
    | otherwise
    |   a.rot.bits = (a.bits + 6) MOD 8

    The arc "Onext(a)" is given by "arc(a.edge.enext[rbit], a.edge.bnext[rbit])",
    provided that the "flip" bit of "a" is 0. Otherwise "Onext(a)" is computed by 
    the formula: "Onext(a) = Flip(Rot(Onext(Rot(Flip(a)))))". *)
    
PROCEDURE OctInit(e: Edge): Edge =    
  BEGIN 
    e.marks[0] := FALSE;
    e.marks[1] := FALSE;
            
    FOR i := 0 TO 3 DO e.enext[i] := e END;
    e.bnext[0] := 0;
    e.bnext[1] := 6;
    e.bnext[2] := 4;
    e.bnext[3] := 2;
    RETURN e
  END OctInit;
  
PROCEDURE MakeEdge (): Arc =
  BEGIN
    WITH e = NEW(Edge).init() DO RETURN Arc{edge := e, bits := 0} END;
  END MakeEdge;
  
PROCEDURE Degree(a: Arc): CARDINAL =
  VAR n: CARDINAL := 0; s: Arc := a;
  BEGIN
    REPEAT INC(n); s := Onext(s) UNTIL s = a;
    RETURN n
  END Degree;

PROCEDURE SenseBit(a: Arc): SBit =
  BEGIN RETURN Word.RightShift(a.bits, 2) END SenseBit;

PROCEDURE FlipBit (a: Arc): FBit =
  BEGIN RETURN Word.And(a.bits, 1) END FlipBit;
  
PROCEDURE DualBit (a: Arc): DBit =
  BEGIN RETURN Word.And(Word.RightShift(a.bits, 1), 1) END DualBit;

PROCEDURE RotBits (a: Arc): RBits =
  BEGIN RETURN Word.And( Word.RightShift(a.bits, 1), 3) END RotBits;
  
PROCEDURE GetArcNum (a: Arc): ArcNum =
  BEGIN RETURN Word.LeftShift(a.edge.num, 3) + a.bits END GetArcNum;

PROCEDURE Rot(a: Arc): Arc =
  BEGIN
    RETURN Arc{edge := a.edge, 
               bits := Word.And( a.bits + 2 + Word.LeftShift( 
               Word.And(a.bits, 1), 2), 7) };
  END Rot;

PROCEDURE Flip(a: Arc): Arc =
  BEGIN
    RETURN Arc{edge := a.edge, bits := Word.Xor( a.bits, 1 )};
  END Flip;
  
PROCEDURE Sym(a: Arc): Arc =
  BEGIN
    RETURN Arc{edge := a.edge, bits := Word.And( (a.bits + 4), 7 )};
  END Sym;
  
PROCEDURE Tor(a: Arc): Arc =
  BEGIN
    RETURN Arc{edge := a.edge,
               bits := Word.And( a.bits + 6 + Word.LeftShift( 
                       Word.And(a.bits, 1), 2), 7) };    
  END Tor;

PROCEDURE Onext(a: Arc): Arc =
  BEGIN
    WITH s = Word.And( Word.RightShift( (a.bits + 1), 1 ), 3 ) DO     
      WITH 
      r = a.edge.bnext[s] 
      DO
        a.edge := a.edge.enext[s];
        IF ( Word.And( a.bits, 1 ) = 1) THEN 
          a.bits := Word.Xor( Word.And( r + 2 + Word.LeftShift( 
                    Word.And(r, 1), 2), 7), 1) ;
        ELSE
          a.bits := r;
        END;
      END
    END;
    RETURN a;
  END Onext;

PROCEDURE Oprev(a: Arc): Arc =
  BEGIN RETURN Rot(Onext(Rot(a))) END Oprev;
  
PROCEDURE Dnext(a: Arc): Arc =
  BEGIN RETURN Sym(Onext(Sym(a)))  END Dnext;

PROCEDURE Dprev(a: Arc): Arc =
  BEGIN RETURN Tor(Onext(Tor(a)))  END Dprev;
  
PROCEDURE Lnext(a: Arc): Arc =
  BEGIN RETURN Rot(Onext(Tor(a)))  END Lnext;
  
PROCEDURE Lprev(a: Arc): Arc =
  BEGIN RETURN Sym(Onext(a))  END Lprev;
  
PROCEDURE Rnext(a: Arc): Arc =
  BEGIN RETURN Tor(Onext(Rot(a)))  END Rnext;

PROCEDURE Rprev(a: Arc): Arc =
  BEGIN RETURN Onext(Sym(a))  END Rprev;

<*UNUSED*>
PROCEDURE Debug() =
  BEGIN 
  END Debug;

PROCEDURE Splice (a, b: Arc) =
  BEGIN
    <* ASSERT b # Flip(Onext(a)) *>
    IF a # b THEN
      WITH 
        ta = Onext(a), tb = Onext(b),
        c = Rot(ta), d = Rot(tb), 
        tc = Onext(c), td = Onext(d)
      DO
        IF ( Word.And(a.bits, 1) = 0 ) THEN 
          WITH
            rba = RotBits(a)
          DO
            a.edge.enext[rba] := tb.edge;
            a.edge.bnext[rba] := tb.bits;
          END
        ELSE
          WITH
            ra = Rot(Flip(a)), fd = Flip(d), rbra = RotBits(ra)
          DO
            ra.edge.enext[rbra] := fd.edge;
            ra.edge.bnext[rbra] := fd.bits;
          END
        END;

        IF ( Word.And(b.bits, 1) = 0 ) THEN 
          WITH
            rbb = RotBits(b)
          DO
            b.edge.enext[rbb] := ta.edge;
            b.edge.bnext[rbb] := ta.bits;          
          END
        ELSE
          WITH
            rb = Rot(Flip(b)), fc = Flip(c), rbrb = RotBits(rb)            
          DO
            rb.edge.enext[rbrb] := fc.edge;
            rb.edge.bnext[rbrb] := fc.bits;
          END
        END;

        IF ( Word.And(c.bits, 1) = 0 ) THEN 
          WITH
            rbc = RotBits(c) 
          DO
            c.edge.enext[rbc] := td.edge;
            c.edge.bnext[rbc] := td.bits;
          END
        ELSE
          WITH
            rc = Rot(Flip(c)), fb = Flip(b), rbrc = RotBits(rc)
          DO
            rc.edge.enext[rbrc] := fb.edge;
            rc.edge.bnext[rbrc] := fb.bits;
          END          
        END;

        IF ( Word.And(d.bits, 1) = 0 ) THEN
        WITH
          rbd = RotBits(d)
        DO
          d.edge.enext[rbd] := tc.edge;
          d.edge.bnext[rbd] := tc.bits;          
        END
        ELSE
          WITH
            rd = Rot(Flip(d)), fa = Flip(a), rbrd = RotBits(rd) 
          DO
            rd.edge.enext[rbrd] := fa.edge;
            rd.edge.bnext[rbrd] := fa.bits;
          END
        END;
      END;
    END;  
  END Splice;

PROCEDURE EnumEdges(a: Arc; visit: VisitProc; edges: BOOLEAN := FALSE)=

  CONST IniStackSize = 1024;

  VAR estack := NEW(REF ARRAY OF Edge, IniStackSize);
      bstack := NEW(REF ARRAY OF FRBits, IniStackSize);
      top: CARDINAL;

  PROCEDURE DoubleStack() =
    BEGIN
      WITH 
        sz = NUMBER(estack^), szNew = 2*sz,
        estackNew = NEW(REF ARRAY OF Edge, szNew),
        bstackNew = NEW(REF ARRAY OF FRBits, szNew)
      DO
        SUBARRAY(estackNew^, 0, sz) := estack^; estack := estackNew;
        SUBARRAY(bstackNew^, 0, sz) := bstack^; bstack := bstackNew;
      END
    END DoubleStack;

  PROCEDURE VisitandMark (c: Arc)=
    (* If edge(c) is unmarked: visit, mark, and stack it. *)
    BEGIN
      IF NOT c.edge.marks[FlipBit(c)] THEN 
        visit(c);
        IF NOT edges THEN visit(Sym(c)) END;
        c.edge.marks[FlipBit(c)] := TRUE; 
        IF top >= NUMBER(estack^) THEN DoubleStack() END;
        estack[top] := c.edge; bstack[top] := c.bits;
        top := top + 1;
      END;
    END VisitandMark;

  VAR seen: CARDINAL; (* # of quads whose childeren were looked at *)
  BEGIN
    <* ASSERT NOT a.edge.marks[FlipBit(a)] *>
    top := 0;
    seen := 0;
    TRY
      VisitandMark (a);
      WHILE seen < top DO 
        WITH b = Arc{edge := estack[seen], bits := bstack[seen]} DO
          VisitandMark(Onext(b));
          VisitandMark(Onext(Sym(b)))
        END;
        seen := seen + 1
      END;
    FINALLY
      (* Erase all marks *)
      WHILE top > 0 DO
        top := top - 1;
        WITH b = Arc{edge := estack[top], bits := bstack[top]} DO
          b.edge.marks[FlipBit(b)] := FALSE;
        END;
      END;
    END
  END EnumEdges;

PROCEDURE EnumVertices (a: Arc; visit: VisitProc)=
  CONST IniStackSize = 1024;

  VAR estack := NEW(REF ARRAY OF Edge, IniStackSize);
      bstack := NEW(REF ARRAY OF FRBits, IniStackSize);
      top: CARDINAL;

  PROCEDURE DoubleStack() =
    BEGIN
      WITH 
        sz = NUMBER(estack^), szNew = 2*sz,
        estackNew = NEW(REF ARRAY OF Edge, szNew),
        bstackNew = NEW(REF ARRAY OF FRBits, szNew)
      DO
        SUBARRAY(estackNew^, 0, sz) := estack^; estack := estackNew;
        SUBARRAY(bstackNew^, 0, sz) := bstack^; bstack := bstackNew;
      END
    END DoubleStack;

  PROCEDURE VisitandMark (c: Arc)=
    (* If org(c) is unmarked: visit, mark, and stack it *)
    VAR cn: Arc;
    BEGIN
      IF NOT c.edge.marks[SenseBit(c)] THEN
        visit(c);
        cn := c;
        REPEAT
          cn.edge.marks[SenseBit(cn)] := TRUE;
          cn := Onext(cn);
        UNTIL (cn = c);
        IF top >= NUMBER(estack^) THEN DoubleStack() END;
        estack[top] := c.edge; bstack[top] := c.bits;
        top := top + 1;
      END;
    END VisitandMark;

  VAR seen: CARDINAL; (* # of quads whose childeren were looked at *)
  BEGIN
    <* ASSERT NOT a.edge.marks[0] *>
    top := 0;
    seen := 0;
    TRY
      VisitandMark (a);
      WHILE seen < top DO
        WITH b = Arc{edge := estack[seen], bits := bstack[seen]} DO
          VAR atn: Arc := b;
          BEGIN
            REPEAT
              VisitandMark(Sym(atn));
              atn := Onext(atn);
            UNTIL (atn = b);
          END
        END;
        seen := seen + 1;
      END; 
    FINALLY
      (* Erase all marks *)
      WHILE top > 0 DO
        top := top - 1;
        WITH b = Arc{edge := estack[top], bits := bstack[top]} DO
          VAR atn := b;
          BEGIN
            REPEAT
              atn.edge.marks[SenseBit(atn)] := FALSE;
              atn:= Onext(atn);
            UNTIL (atn = b);
          END
        END
      END;
    END;
  END EnumVertices;
  
PROCEDURE NumberEdges(READONLY a: ARRAY OF Arc): REF ARRAY OF Arc =
  CONST IniStackSize = 1024;

  VAR stack := NEW(REF ARRAY OF Arc, IniStackSize);
      nstack: CARDINAL := 0;

  PROCEDURE DoubleStack() =
    BEGIN
      WITH 
        sz = NUMBER(stack^), szNew = 2*sz,
        stackNew = NEW(REF ARRAY OF Arc, szNew)
      DO
        SUBARRAY(stackNew^, 0, sz) := stack^; stack := stackNew;
      END
    END DoubleStack;

      
  (* An edge "e" is "marked" if stack[e.num].edge = e. *)

  PROCEDURE VisitAndMark (t: Arc) =
    (* If t is unmarked: visit, mark, and stack it. *)
    BEGIN
      WITH tnum = t.edge.num DO
        IF tnum < nstack AND stack[tnum].edge = t.edge THEN
          (* Edge is already marked *)
        ELSE
          IF nstack = NUMBER(stack^) THEN DoubleStack() END;
          tnum := nstack;
          stack[nstack] := t;
          INC(nstack);
        END
      END
    END VisitAndMark;
  
  VAR seen: CARDINAL := 0; (* # of edges whose childeren were looked at *)
  BEGIN
    nstack := 0;
    seen := 0;
    FOR i := 0 TO LAST(a) DO
      VisitAndMark (a[i])
    END;
    WHILE seen < nstack DO
      WITH s = stack[seen] DO
        VisitAndMark(Onext(s));
        VisitAndMark(Onext(Sym(s)))
      END;
      seen := seen + 1
    END;
    WITH
      r = NEW(REF ARRAY OF Arc, nstack)
    DO
      r^ := SUBARRAY(stack^, 0, nstack);
      RETURN r
    END;
  END NumberEdges;
  
PROCEDURE SetOnext(a, b: Arc) =
  BEGIN
    IF Onext(a) # b THEN Splice(a, Oprev(b)) END
  END SetOnext;

PROCEDURE PrintArc (wr: Wr.T; a: Arc; eWidth: CARDINAL := 1) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText (wr, 
      Fmt.Pad(Fmt.Int(a.edge.num), eWidth) & ":" &
      Fmt.Int(RotBits(a)) & ":" &
      Fmt.Int(FlipBit(a))
    )
  END PrintArc;
  
PROCEDURE ReadArc(rd: Rd.T; READONLY map: ARRAY OF Edge): Arc =
  VAR f, n, r: CARDINAL;
  <* FATAL Rd.Failure, Rd.EndOfFile, Thread.Alerted, FloatMode.Trap, Lex.Error *>
  BEGIN
    Lex.Skip(rd);
    n := Lex.Int(rd);
    <* ASSERT Rd.GetChar(rd) = ':' *>
    r := Lex.Int(rd);
    <* ASSERT Rd.GetChar(rd) = ':' *>
    f := Lex.Int(rd);
    <* ASSERT n < NUMBER(map) *>
    <* ASSERT r < 4 AND f < 2 *>
    RETURN Arc{edge := map[n], bits := 2 * r + f}
  END ReadArc;

PROCEDURE PrintEdge(wr: Wr.T; e: Edge; eWidth: CARDINAL := 1) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    FOR i := 0 TO 3 DO 
      IF i > 0 THEN Wr.PutChar(wr, ' ') END;
      PrintArc(wr, Arc{edge := e.enext[i], bits := e.bnext[i]}, eWidth);
    END
  END PrintEdge;
  
PROCEDURE ReadEdge(rd: Rd.T; e: Edge; READONLY map: ARRAY OF Edge) =
  BEGIN
    FOR i := 0 TO 3 DO 
      SetOnext(Arc{edge := e, bits := 2 * i}, ReadArc(rd, map))
    END
  END ReadEdge;
  
BEGIN
END Oct.
 
