(* Last edited on 1999-06-05 20:50:59 by stolfi *)

MODULE DAG;

(* See the copyright and disclaimer note at the end of this file. *)

IMPORT Rd, Wr, Thread, Text;
IMPORT FileFmt, OldFileFmt, FGet, NGet;
IMPORT BinGraph, BinGraphF, DAGF;  (* Reveals internal structure of a BinGraph.T *)
FROM Basics IMPORT NAT;
FROM Basics IMPORT Skip, Abort;

REVEAL
  T = DAGF.Private BRANDED OBJECT
    OVERRIDES
      (* Inherted from "DAG.Public": *)
      Class     := Class;
      IsSink    := IsSink;
      Last      := Last;
      Rest      := Rest;
      Load      := Load;
      Dump      := Dump;
      First     := First;
      OutDeg    := OutDeg;
      EnumPaths := EnumPaths;

      (* Inherited from "DAGF.Private" *)
      Discard   := Discard;
    END;

TYPE 
  Entry = BinGraphF.Entry;

PROCEDURE Class(dag: T; s: Node): Symbol =
  BEGIN
    <* ASSERT s < dag.nNodes *>
    WITH es = dag.e^[s] DO
      RETURN es.lSym
    END;
  END Class;

PROCEDURE IsSink(dag: T; s: Node): BOOLEAN =
  BEGIN
    <* ASSERT s < dag.nNodes *>
    WITH es = dag.e[s] DO
      RETURN es.rNode # s
    END;
  END IsSink;

PROCEDURE Last(dag: T; s: Node): ArcData =
  BEGIN
    <* ASSERT s < dag.nNodes *>
    WITH es = dag.e[s] DO
      <* ASSERT es.lNode < s *>
      <* ASSERT es.rNode < s *>
      RETURN ArcData{label := es.rSym, dest := es.rNode}
    END
  END Last;

PROCEDURE Rest(dag: T; s: Node): Node =
  BEGIN
    <* ASSERT s < dag.nNodes *>
    WITH es = dag.e[s] DO
      <* ASSERT es.lNode < s *>
      <* ASSERT es.rNode < s *>
      RETURN es.lNode
    END;
  END Rest;

PROCEDURE OutDeg(dag: T; s: Node): NAT =
  VAR t: Node := s; n: NAT := 0;
  BEGIN
    WITH e = dag.e^ DO
      LOOP
        WITH et = e[t] DO
          IF et.rNode = t THEN RETURN n END;
          INC(n); t := et.lNode
        END
      END
    END;
  END OutDeg;

PROCEDURE First(dag: T; s: Node): ArcData =
  VAR t: Node := s; p: Node := NoNode;
  BEGIN
    WITH e = dag.e^ DO
      LOOP
        WITH et = e[t] DO
          IF et.rNode = t THEN EXIT END;
          p := t; t := et.lNode
        END
      END;
      <* ASSERT p # NoNode *>
      WITH ep = e[p] DO
        RETURN ArcData{label := ep.rSym, dest := ep.rNode}
      END
    END;
  END First;

PROCEDURE CheckInvariants(dag: T; lo, hi: Node) =
  (*
    Checks representation invariants for all nodes 
    "s" in "[lo .. hi-1]". *)
  BEGIN
    IF dag.e = NIL OR lo > dag.nNodes THEN RETURN END;
    WITH e = dag.e^ DO
      FOR s := lo TO MIN(hi, dag.nNodes)-1 DO
        WITH es = e[s] DO
          IF es.lNode = s THEN
            (* Must be a sink: *)
            <* ASSERT es.rNode = s *>
            <* ASSERT es.rNode = 0 *>
          ELSE
            (* Must be topologically sorted,and class-homogeneous: *)
            <* ASSERT es.rNode < s *>
            <* ASSERT es.lNode < s *>
            <* ASSERT e[es.lNode].lSym = es.lSym *>
          END
        END
      END
    END
  END CheckInvariants;

PROCEDURE Discard(dag: T; s: Node): T =
  BEGIN
    <* ASSERT s < dag.nNodes *>
    dag.nNodes := s;
    INC(dag.epoch);
    RETURN dag
  END Discard;

PROCEDURE EnumPaths(
    dag: T;
    s: Node;
    enter: NodeAction := NIL;
    push, pop: ArcAction := NIL;
    exit: NodeAction := NIL;
  ) RAISES {Abort} =

  VAR len: NAT := 0;

  BEGIN
  
    <* ASSERT s < dag.nNodes *>
    
    WITH e = dag.e^ DO

      PROCEDURE DoEnumPaths(t: Node) RAISES {Abort} =
        (*
          Does EnumPaths starting at "t", a generic sucessor of "s". *)

        VAR i: NAT := 0;

        PROCEDURE EnumRest(r: Node) RAISES {Skip, Abort} =
          (*
            Enumerates a prefix "r" OF node "t". *)

          BEGIN
            WITH er = e[r] DO 
              IF er.rNode # r THEN EnumRest(er.lNode) END;
              TRY
                IF push # NIL THEN push(len, t, er.lSym, i, er.rSym, er.rNode) END;
                INC(len);
                DoEnumPaths(er.rNode);
                DEC(len);
              EXCEPT
                Skip => (* Ok *)
              END;
              IF pop # NIL THEN pop(len, t, er.lSym, i, er.rSym, er.rNode) END;
              INC(i);
            END
          END EnumRest;

        BEGIN
          WITH et = e[t] DO
            TRY 
              IF enter # NIL THEN enter(len, t, et.lSym) END;
              EnumRest(t);
            EXCEPT Skip => (* Ok *) END;
            IF exit # NIL THEN
              <* FATAL Skip *> BEGIN exit(len, t, et.lSym) END
            END
          END;
        END DoEnumPaths;

      BEGIN
        <* ASSERT s < dag.nNodes *>
        DoEnumPaths(s)
      END

    END
  END EnumPaths;

CONST
  FileType    = "DAG.T";
  FileVersion = "97-01-21";

  OldFileType     = "DAG.Dump";
  OldFileVersion  = "91-12-16";

PROCEDURE Dump(dag: T; wr: Wr.T) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    FileFmt.WriteHeader(wr, FileType, FileVersion);
    NARROW(dag, BinGraph.T).Dump(wr);
    FileFmt.WriteFooter(wr, FileType);
    Wr.Flush(wr);
  END Dump;

PROCEDURE Load(rd: Rd.T): T =
  <* FATAL Rd.Failure, Rd.EndOfFile, Thread.Alerted *>
  BEGIN
    WITH dag = NEW(T) DO
      (* Check file version *)
      WITH hdr = Rd.GetLine(rd) DO
        IF Text.Equal(hdr, FileFmt.MakeHeader(FileType, FileVersion)) THEN
          dag := NARROW(dag, BinGraph.T).LoadBody(rd);
          FileFmt.ReadFooter(rd, FileType)
        ELSIF Text.Equal(hdr, OldFileFmt.MakeHeader(OldFileType, OldFileVersion)) THEN
          dag.comment := FileFmt.ReadComment(rd, '|');
          dag := NARROW(dag, BinGraph.T).LoadOldDAGBody(rd);
          OldFileFmt.ReadFooter(rd, OldFileType)
        ELSE
          <* ASSERT FALSE *>
        END;
      END;
      CheckInvariants(dag, 0, dag.nNodes);
      RETURN dag
    END
  END Load;

BEGIN
END DAG.

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
