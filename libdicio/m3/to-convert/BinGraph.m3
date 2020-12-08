(* Last edited on 1999-06-05 20:50:50 by stolfi *)

MODULE BinGraph;

(* See the copyright and disclaimer note at the end of this file. *)

IMPORT Rd, Wr, Thread;
IMPORT FileFmt, FGet, FPut, NGet, NPut;
IMPORT BinGraphF;
FROM Basics IMPORT NAT;
FROM Basics IMPORT Full, Skip, Abort;
FROM BinGraphF IMPORT Entry, Entries; 

REVEAL
  T = BinGraphF.Private BRANDED OBJECT
    OVERRIDES
      (* Inherited from "BinGraph.Public" *)
      NNodes    := NNodes;
      Get       := Get;
      EnumPaths := EnumPaths;
      Dump      := Dump;
      Load      := Load;

      (* Inherited from "BinGraphF.Private" *)
      NAlloc    := NAlloc;
      Alloc     := Alloc;
      Set       := Set;
      Add       := Add;
      Clear     := Clear;
      Copy      := Copy;
      Crunch    := Crunch;
    END;

PROCEDURE NNodes(g: T): Node =
  BEGIN
    RETURN g.nNodes
  END NNodes;

PROCEDURE Get(g: T; s: Node): NodeData =
  BEGIN
    <* ASSERT s < g.nNodes *>
    WITH es = g.e[s] DO
      RETURN NodeData{
        lSym := es.lSym, lNode := es.lNode, 
        rSym := es.rSym, rNode := es.rNode
      }
    END
  END Get;

PROCEDURE NAlloc(g: T): NAT =
  BEGIN
    IF g.e = NIL THEN RETURN 0 ELSE RETURN NUMBER(g.e^) END;
  END NAlloc;

PROCEDURE Alloc(g: T; maxNodes: NAT): T =
  BEGIN
    <* ASSERT maxNodes >= g.nNodes *>
    (* See if we already have the right amount of storage: *)
    WITH nAlloc = g.NAlloc() DO
      IF nAlloc >= maxNodes AND nAlloc - maxNodes < maxNodes THEN RETURN g END;
      g.e := NIL;
      IF maxNodes = 0 THEN RETURN g END;
    END;
    
    (* Needs to (re)allocate: *)
    WITH
      e = g.e,
      eNew = NEW(REF Entries, maxNodes)
    DO
      IF g.nNodes > 0 THEN
        SUBARRAY(eNew^, 0, g.nNodes) := e^
      END;
      e := eNew
    END;
    RETURN g
  END Alloc;

PROCEDURE Set(g: T; s: Node; lSym, rSym: Symbol; lNode, rNode: Node) =
  BEGIN
    <* ASSERT s < g.nNodes *>
    <* ASSERT lNode < g.nNodes *>
    <* ASSERT rNode < g.nNodes *>
    WITH e = g.e^ DO
      WITH 
        es = e[s],
        en = Entry{
          es.lSym := lSym,
          es.lNode := lNode,
          es.rSym := rSym,
          es.rNode := rNode
        }
      DO
        IF en # es THEN 
          es := en;
          INC(g.epoch)
        END
      END
    END
  END Set;
  
PROCEDURE Add(g: T; lSym, rSym: Symbol; lNode, rNode: Node): Node RAISES {Full} =
  BEGIN
    <* ASSERT s <= g.nNodes *>
    IF g.e = NIL OR s > LAST(g.e^) THEN RAISE Full END;
    WITH 
      e = g.e^,
      s = g.nNodes + 0
    DO
      INC(g.nNodes);
      <* ASSERT lNode < g.nNodes *>
      <* ASSERT rNode < g.nNodes *>
      WITH es = e[s],
        en = Entry{
          es.lSym := lSym,
          es.lNode := lNode,
          es.rSym := rSym,
          es.rNode := rNode
        }
      DO
        es := en
      END;
      RETURN s
    END
  END Add;
  
PROCEDURE Clear(g: T): T =
  BEGIN
    IF g.nNodes > 0 THEN
      g.nNodes := 0;
      g.e := NIL;
      INC(g.epoch)
    END;
    RETURN g
  END Clear;

PROCEDURE Copy(g: T; from: T; s: Node; VAR map: ARRAY OF Node): Node RAISES {Full} =

  PROCEDURE DoCopy(u: Node): Node RAISES {Full} =
    (*
      Does Copy starting at "u", a generic sucessor of "s": copies
      "u.lNode", then copies "u.rNode", and finally creates the copy of "u". *)

    BEGIN
      IF map[u] = NoNode THEN
        map[u] := g.nNodes;
        WITH 
          eu = from.e[u],
          lNode = DoCopy(eu.lNode),
          rNode = DoCopy(eu.rNode)
        DO
          g.Set(map[u], eu.lSym, eu.rSym, lNode, rNode);
        END
      END;
      RETURN map[u]
    END DoCopy;

  BEGIN
    RETURN DoCopy(s)
  END Copy;

PROCEDURE Crunch(
    g: T; 
    READONLY root: ARRAY OF Node; 
    VAR (*OUT*) map: ARRAY OF Node;
  ): T =
  VAR nOld: NAT;
      nNew: NAT;
  BEGIN
    (*
      Crunch is done in in four passes:

      First pass: unmark all nodes ("map[s] = NoNode").

      Second pass: recursively mark ("map[s] = s") all 
      nodes reachable from "root".

      Third pass: scan all entries from 1 up, and move each marked
      node "s" to the lowest unused position, updating its "l" and
      "r" pointers, and saving that position in "map[s]".  Also
      update the root pointers.

      Fourth pass: rebuild the hash table for the compacted nodes.
    *)

    IF g.e = NIL THEN g.e := NEW(REF Entries, 0) END;
    WITH e = g.e^ DO

      (* Unmark all nodes *)
      FOR i := 0 TO g.nNodes-1 DO map[i] := NoNode END;

      (* Mark reachable nodes, and set "nOld > t" for all copied nodes "t" : *)
      (* (convention: a node "t" is marked iff "map[t]=t".) *)
      
      nOld := 0;
      
      PROCEDURE DoMark(t: Node) =
        BEGIN
          WHILE map[t] = NoNode DO
            IF t >= nOld THEN nOld := t+1 END;
            map[t] := t;
            WITH lNode = e[t].lNode, rNode = e[t].rNode DO
              IF lNode < rNode THEN
                DoMark(lNode); t := rNode
              ELSE
                DoMark(rNode); t := lNode
              END
            END
          END
        END DoMark;
        
      BEGIN
        FOR k := 0 TO LAST(root) DO DoMark(root[k]) END
      END;
      
      (* Compress nodes virtually, updating "map": *)
      nNew := 0;
      FOR u := 0 TO nOld-1 DO
        IF map[u] # NoNode THEN
          map[u] := nNew; INC(nNew)
        END
      END;
      
      (* Now copy the node data, updating links: *)
      FOR u := 0 TO nOld-1 DO
        IF map[u] # NoNode THEN
          WITH
            eu = e[u], 
            lNode = map[eu.lNode], 
            rNode = map[eu.rNode],
            v = map[u], 
            ev = e[v]
          DO
            <* ASSERT lNode # NoNode *>
            <* ASSERT rNode # NoNode *>
            ev.lSym := eu.lSym;
            ev.rSym := eu.rSym;
            ev.lNode := lNode;
            ev.rNode := rNode
          END;
        END
      END;
      IF nNew < g.nNodes THEN
        INC(g.epoch);
        g.nNodes := nNew
      END
    END;
    RETURN g
  END Crunch;

PROCEDURE EnumPaths(
    g: T;
    s: Node;
    enter: NodeAction := NIL;
    push, pop: LinkAction := NIL;
    exit: NodeAction := NIL;
  ) RAISES {Abort} =

  VAR len: NAT := 0;

  PROCEDURE DoEnumPaths(o: Node) RAISES {Abort} =
    (*
      Does "EnumPaths" starting at "o", a generic sucessor of "s". *)
    VAR t: Node; x: Symbol;
    BEGIN
      WITH eo = g.e[o] DO
        TRY
          IF enter # NIL THEN enter(len, o, eo.lSym, eo.rSym) END;
          FOR b := FALSE TO TRUE DO
            IF NOT b THEN 
              t := eo.lNode; x := eo.lSym
            ELSE
              t := eo.rNode; x := eo.lSym
            END;
            TRY
              IF push # NIL THEN push(len, o, b, x, t) END;
              INC(len);
              DoEnumPaths(t);
              DEC(len)
            EXCEPT 
              Skip => (*OK*)
            END;
            IF pop # NIL THEN pop(len, o, b, x, t) END;
          END;
        EXCEPT
          Skip => (*OK*)
        END;
        IF exit # NIL THEN 
          <*FATAL Skip *> BEGIN exit(len, o, eo.lSym, eo.rSym) END
        END
      END
    END DoEnumPaths;

  BEGIN
    <* ASSERT s < g.nNodes *>
    <* ASSERT enter # NIL OR push # NIL *>
    DoEnumPaths(s)
  END EnumPaths;

CONST
  FileType = "BinGraph.T";
  FileVersion = "97-01-21";

  OldFileType     = "DAG.Dump";
  OldFileVersion  = "91-12-16";

PROCEDURE Dump(g: T; wr: Wr.T) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    FileFmt.WriteHeader(wr, FileType, FileVersion);
    FileFmt.WriteComment(wr, g.comment, '|');
    DumpBody(g, wr);
    FileFmt.WriteFooter(wr, FileType);
    Wr.Flush(wr);
  END Dump;

PROCEDURE DumpBody(g: T; wr: Wr.T) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    NPut.Int(wr, "nNodes", g.nNodes);
    WITH e = g.e^ DO
      FOR i := 0 TO g.nNodes-1 DO
        WITH ei = e[i] DO
          FPut.Int(wr, i);        FPut.Space(wr, 1);
          FPut.Int(wr, ei.lSym);  FPut.Space(wr, 1);
          FPut.Int(wr, ei.lNode); FPut.Space(wr, 1);
          FPut.Int(wr, ei.rSym);  FPut.Space(wr, 1);
          FPut.Int(wr, ei.rNode); FPut.EOL(wr);
        END;
      END;
    END;
    Wr.Flush(wr);
  END DumpBody;

PROCEDURE Load(g: T; rd: Rd.T): T =
  <* FATAL Rd.Failure, Rd.EndOfFile, Thread.Alerted *>
  BEGIN
    WITH hdr = Rd.GetLine(rd) DO
      g.comment := FileFmt.ReadComment(rd, '|');
      IF Text.Equal(hdr, FileFmt.MakeHeader(FileType, FileVersion)) THEN
        (* Official format. *)
        LoadBody(g, rd);
        FileFmt.ReadFooter(rd, FileType)
      ELSIF Text.Equal(hdr, OldFileFmt.MakeHeader(OldFileType, OldFileVersion)) THEN
        (* Old "Dag.Dump" format. *)
        LoadOldDAGBody(g, rd);
        OldFileFmt.ReadFooter(rd, OldFileType);
      END;
    END;
    RETURN g
  END Load;

PROCEDURE LoadBody(g: T; rd: Rd.T) =
  BEGIN
    WITH nNodes = NGet.Int(rd, "nNodes") DO
      (* Discard old contents and reallocate if necessary: *)
      IF g.nNodes > 0 THEN
        g.nNodes := 0;
        INC(g.epoch)
      END;
      EVAL g.Alloc(nNodes);
      (* Read nodes: *)
      WITH e = g.e^ DO
        FOR i := 0 TO nNodes-1 DO
          WITH
            s = FGet.Int(rd),
            lSym  = FGet.Int(rd),
            lNode = FGet.Int(rd),
            rSym  = FGet.Int(rd),
            rNode = FGet.Int(rd)
          DO
            <* ASSERT s = i *>
            <* ASSERT lNode < nNodes *>
            <* ASSERT rNode < nNodes *>
            e[s] := Entry{
              lSym  := lSym,
              lNode := lNode,
              rSym  := rSym,
              rNode := rNode
            };
          END;
          FGet.EOL(rd);
        END;
      END;
      g.nNodes := nNodes;
    END
  END LoadBody;

PROCEDURE LoadOldDAGBody(g: T; rd: Rd.T) =
  BEGIN
    WITH nNodes = NGet.Int(rd, "max state") + 1 DO
      (* Discard old contents and reallocate if necessary: *)
      IF g.nNodes > 0 THEN
        g.nNodes := 0;
        INC(g.epoch)
      END;
      EVAL g.Alloc(nNodes);
      WITH e = g.e^ DO
        (* Create the "NullState": *)
        e[0] := Entry{lSym := 0, lNode := 0, rSym := 0, rNode := 0};
        (* Read proper states: *)
        FOR i := 1 TO nNodes-1 DO
          WITH
            s = FGet.Int(rd),
            rSym = FGet.Int(rd),
            lSym = FGet.Int(rd),
            rNode = FGet.Int(rd),
            lNode = FGet.Int(rd)
          DO
            <* ASSERT s = i *>
            <* ASSERT lNode < nNodes *>
            <* ASSERT rNode < nNodes *>
            e[s] := Entry{
              lSym  := lSym,
              lNode := lNode,
              rSym  := rSym,
              rNode := rNode
            }
          END;
          FGet.EOL(rd)
        END
      END;
      g.nNodes := nNodes;
    END
  END LoadOldDAGBody;

BEGIN
END BinGraph.

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
