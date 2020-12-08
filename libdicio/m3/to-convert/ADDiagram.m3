(****************************************************************************)
(* (C) Copyright 1992 Universidade Estadual de Campinas (UNICAMP)           *)
(*                    Campinas, SP, Brazil                                  *)
(*                                                                          *)
(* Authors:                                                                 *)
(*                                                                          *)
(*   Tomasz Kowaltowski  - CS Dept, UNICAMP <tomasz@dcc.unicamp.ansp.br>    *)
(*   Claudio L. Lucchesi - CS Dept, UNICAMP <lucchesi@dcc.unicamp.ansp.br>  *)
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

MODULE ADDiagram;

FROM ADTypes IMPORT 
  NAT, BOOL, IntPoint, IntBox, IntRange, Counts, 
  RealPoint, ArcControlPoints, Nowhere, ClipIntPoint;
FROM ADMath IMPORT Sin, Cos, Dir, Tan, Sqrt, Sqr;
FROM ADReport IMPORT 
  ReportUnplacement, ReportPlacement,
  ReportUnhappiness, ReportArcArcCollision,
  ReportEnergyComputation, ReportArcStrain;
  
TYPE
  BoolMap = REF ARRAY OF ARRAY OF BOOL;        (* Marks grid cells *)
  ArcNodeMap = REF ARRAY OF ARRAY OF Fragment; (* Maps cell indices to diagram fragments *)
 
REVEAL
  T = Public BRANDED OBJECT
      who: ArcNodeMap;       (* Which nodes/arcs are using cell [x,y]. *)
      seen: BoolMap;         (* Which cells have already been enumerated *)
      busy: BOOL;            (* Enumeration in progress. *)
    OVERRIDES
      AddNode := AddNode;
      EnumNodes := EnumNodes;
      EnumNodeCells := EnumNodeCells;
      AddArc := AddArc;
      EnumArcs := EnumArcs;
      EnumArcCells := EnumArcCells;
      MoveNode := MoveNode;
      UnplaceNode := UnplaceNode;
      PlaceNode := PlaceNode;
      MakeNodeUnhappy := MakeNodeUnhappy;
    END;
    
  Element = Fragment BRANDED OBJECT
    OVERRIDES
      AddElement := ElementAddElement;
      SubElement := ElementSubElement;
      Includes := ElementIncludes;
    END;

  Node = NodePublic BRANDED OBJECT
    OVERRIDES
      EnumElements := NodeEnumElements;
    END;
    
  Arc = ArcPublic BRANDED OBJECT
    OVERRIDES
      EnumElements := ArcEnumElements;
    END;
    
TYPE
  FragmentPair = Fragment BRANDED OBJECT
      a, b: Fragment; 
    OVERRIDES
      EnumElements := PairEnumElements;
      AddElement := PairAddElement;
      SubElement := PairSubElement;
      Includes := PairIncludes;
    END;
    
PROCEDURE NodeEnumElements(
    n: Node; 
    <*UNUSED*> arcAction: ArcAction; 
    nodeAction: NodeAction;
  ) RAISES {Abort} =
  BEGIN
    IF nodeAction # NIL THEN nodeAction(n) END
  END NodeEnumElements;

PROCEDURE ArcEnumElements(
    a: Arc; 
    arcAction: ArcAction; 
    <*UNUSED*> nodeAction: NodeAction;
  ) RAISES {Abort} =
  BEGIN
    IF arcAction # NIL THEN arcAction(a) END
  END ArcEnumElements;

PROCEDURE PairEnumElements(
    p: FragmentPair; 
    arcAction: ArcAction; 
    nodeAction: NodeAction;
  ) RAISES {Abort} =
  BEGIN
    p.a.EnumElements(arcAction, nodeAction);
    p.b.EnumElements(arcAction, nodeAction);
  END PairEnumElements;

PROCEDURE ElementAddElement (e: Element; x: Element): Fragment =
  BEGIN
    IF x = e THEN 
      RETURN e
    ELSE
      RETURN NEW(FragmentPair, a := e, b := x)
    END
  END ElementAddElement;

PROCEDURE PairAddElement (p: FragmentPair; x: Element): Fragment =
  BEGIN
    IF p.a.Includes(x) OR p.b.Includes(x) THEN
      RETURN p
    ELSE
      RETURN NEW(FragmentPair, a := p.b, b := p.a.AddElement(x))
    END
  END PairAddElement;

PROCEDURE ElementSubElement (e: Element; x: Element): Fragment =
  BEGIN
    IF x = e THEN 
      RETURN NIL
    ELSE
      RETURN e
    END
  END ElementSubElement;

PROCEDURE PairSubElement (p: FragmentPair; x: Element): Fragment =
  VAR na, nb: Fragment;
  BEGIN
    nb := p.b.SubElement(x);
    IF nb = p.b THEN na := p.a.SubElement(x) ELSE na := p.a END;
    IF na = p.a AND nb = p.b THEN
      RETURN p
    ELSIF na = NIL THEN
      RETURN nb
    ELSIF nb = NIL THEN
      RETURN na
    ELSE
      RETURN NEW(FragmentPair, a := nb, b := na)
    END
  END PairSubElement;

PROCEDURE ElementIncludes (e: Element; x: Element): BOOL =
  BEGIN
    RETURN x = e
  END ElementIncludes;

PROCEDURE PairIncludes (p: FragmentPair; x: Element): BOOL =
  BEGIN
    RETURN p.b.Includes(x) OR p.a.Includes(x);
  END PairIncludes;

PROCEDURE New(maxLabel: CARDINAL; maxNodes: CARDINAL; gridN: Counts): T =
  BEGIN
    WITH
      t = NEW(T,
        nNodes := 0,
        node := NEW(REF ARRAY OF Node, maxLabel + 1),
        gridN := gridN,
        unplaced := NewQueue(maxNodes),
        unhappy := NewQueue(maxNodes),
        who := NEW(ArcNodeMap, gridN[0], gridN[1]),
        seen := NEW(BoolMap, gridN[0], gridN[1]),
        busy := FALSE
      )
    DO
      (* Clear labe-to-node map *)
      WITH node = t.node^ DO
        FOR i := 0 TO maxLabel DO node[i] := NIL END;
      END;
      (* Clear "who" and "seen" maps: *)
      WITH who = t.who^, seen = t.seen^ DO
        FOR x := 0 TO gridN[0] - 1 DO
          FOR y := 0 TO gridN[1] - 1 DO
            who[x,y] := NIL;
            seen[x, y] := FALSE
          END
        END;
      END;
      RETURN t
    END;
  END New;

PROCEDURE AddNode(
    t: T;
    label: NodeLabel; (* Node label *)
    marked: BOOL; (* TRUE to highlight the node *)
    ri: NAT;      (* Radius of input fan, in grid cells *)
    ro: NAT;      (* Radius of output fan, in grid cells *)
    ht: NAT;      (* Half-height, in grid cells *)
  ): Node =
  BEGIN
    <* ASSERT label <= LAST(t.node^) *>
    WITH 
      si = NodeSideData{rad := ri, arcs := NIL},
      so = NodeSideData{rad := ro, arcs := NIL},
      n = NEW(Node,
        label := label,
        gp := Nowhere,
        range := IntBox{IntRange{0, t.gridN[0]}, IntRange{0, t.gridN[1]}},
        placed := FALSE,
        happy := FALSE,
        marked := marked,
        when := 0,
        ht := ht,
        tilt := 0.0,
        s := ARRAY [0..1] OF NodeSideData{si, so}
      ),
      node = t.node^
    DO
      INC(t.nNodes);
      <* ASSERT node[label] = NIL *>
      node[label] := n;
      PutNodeInQueue(t.unplaced, n);
      RETURN n
    END;
  END AddNode;
  
PROCEDURE EnumNodes(t: T; action: NodeAction) RAISES {Abort} =
  BEGIN
    WITH node = t.node^ DO
      FOR i := 0 TO LAST(node) DO
        IF node[i] # NIL THEN action(node[i]) END
      END
    END
  END EnumNodes;

PROCEDURE EnumNodeCells(t: T; n: Node; action: NodeCellAction) RAISES {Abort} =
  BEGIN
    WITH
      gridN = t.gridN,
      gp = n.gp,
      who = t.who^,
      ri = n.s[0].rad,
      ro = n.s[1].rad,
      ht = n.ht
    DO
      FOR i0 := MAX(0, gp[0] - ri) TO MIN(gridN[0], gp[0] + ro) - 1 DO
        FOR i1 := MAX(0, gp[1] - ht) TO MIN(gridN[1], gp[1] + ht) - 1 DO
          action(n, i0, i1, who[i0, i1])
        END
      END;
    END
  END EnumNodeCells;
  
PROCEDURE NodeBBox(t: T; n: Node): IntBox =
  BEGIN
    WITH
      gridN = t.gridN,
      gp = n.gp,
      ri = n.s[0].rad,
      ro = n.s[1].rad,
      ht = n.ht
    DO
      RETURN IntBox{
        IntRange{gp[0] - ri, gp[0] + ro - 1},
        IntRange{gp[1] - ht, gp[1] + ht - 1}
      }
    END
  END NodeBBox;

PROCEDURE AddArc(
    <*UNUSED*> t: T;
    label: ArcLabel; (* Arc label *)
    org: Node;       (* Origin node *)
    dst: Node;       (* Destination node *)
  ): Arc =
  BEGIN
    <* ASSERT org # NIL *>
    <* ASSERT dst # NIL *>
    WITH
      so = ArcEndData{node := org, next := org.s[1].arcs, dir := 0.0},
      sd = ArcEndData{node := dst, next := dst.s[0].arcs, dir := 0.0},
      a = NEW(Arc,
        label := label,
        s := ARRAY [0..1] OF ArcEndData{so, sd}
      )
    DO
      org.s[1].arcs := a;
      dst.s[0].arcs := a;
      RETURN a
    END
  END AddArc;

PROCEDURE EnumArcs(t: T; action: ArcAction) RAISES {Abort} =
  BEGIN
    WITH node = t.node^ DO
      FOR i := 0 TO LAST(node) DO
        IF node[i] # NIL THEN
          VAR a: Arc := node[i].s[1].arcs;
          BEGIN
            WHILE a # NIL DO
              action(a);
              a := a.s[0].next
            END
          END;
        END;
      END
    END
  END EnumArcs;

PROCEDURE EnumArcCells(t: T; a: Arc; action: ArcCellAction) RAISES {Abort} =
  BEGIN
    WITH
      gridN = t.gridN,
      orgBB = NodeBBox(t, a.s[0].node),
      dstBB = NodeBBox(t, a.s[1].node),
      cp = ComputeArcControlPoints(a)
    DO
      <* ASSERT NOT t.busy *>
      TRY
        t.busy := TRUE;
        ClearBoolMap(t.seen^, SegmentGridBox(cp[0], cp[1], gridN));
        ClearBoolMap(t.seen^, BezierGridBox(cp[1], cp[2], cp[3], cp[4], gridN));
        EnumSegmentCells(t, a, cp[0], cp[1], orgBB, dstBB, action);
        EnumBezierCells(t, a, cp[1], cp[2], cp[3], cp[4], orgBB, dstBB, action)
      FINALLY
        t.busy := FALSE
      END
    END
  END EnumArcCells;
  
PROCEDURE ClearBoolMap(VAR map: ARRAY OF ARRAY OF BOOL; READONLY box: IntBox) =
  BEGIN
    FOR i0 := box[0].lo TO box[0].hi DO
      FOR i1 := box[1].lo TO box[1].hi DO
        map[i0, i1] := FALSE
      END
    END;
  END ClearBoolMap;
  
PROCEDURE EnumSegmentCells(
    t: T; 
    a: Arc; 
    p0, p1: RealPoint; 
    READONLY orgBB, dstBB: IntBox;
    action: ArcCellAction;
  ) RAISES {Abort} =
  
  BEGIN
    WITH
      gridN = t.gridN,
      who = t.who^,
      seen = t.seen^
    DO

      PROCEDURE DoEnum(READONLY p0, p1: RealPoint) RAISES {Abort} =
        BEGIN
          WITH
            box = SegmentGridBox(p0, p1, t.gridN),
            size0 = box[0].hi - box[0].lo + 1,
            size1 = box[1].hi - box[1].lo + 1,
            minSize = MIN(size0, size1),
            maxSize = MAX(size0, size1)
          DO
            IF maxSize <= 2  OR minSize <= 1 THEN 
              FOR i0 := box[0].lo TO box[0].hi DO
                FOR i1 := box[1].lo TO box[1].hi DO
                  IF NOT InBox(i0, i1, orgBB) 
                  AND NOT InBox(i0, i1, dstBB) 
                  AND NOT t.seen[i0, i1] THEN
                    t.seen[i0, i1] := TRUE;
                    action(a, i0, i1, t.who[i0, i1])
                  END;
                END
              END;
              RETURN
            END
          END;
          WITH pm = Interpolate(p0, p1, 0.5) DO
            DoEnum(p0, pm);
            DoEnum(pm, p1)
          END
        END DoEnum;

      BEGIN
        DoEnum(p0, p1);
      END
    END
  END EnumSegmentCells;
  
PROCEDURE EnumBezierCells(
    t: T; 
    a: Arc; 
    p000, p001, p011, p111: RealPoint; 
    READONLY orgBB, dstBB: IntBox;
    action: ArcCellAction;
  ) RAISES {Abort} =
  BEGIN
    WITH
      gridN = t.gridN,
      who = t.who^,
      seen = t.seen^
    DO

      PROCEDURE DoEnum(READONLY p000, p001, p011, p111: RealPoint) RAISES {Abort} =
        BEGIN
          WITH
            box = BezierGridBox(p000, p001, p011, p111, gridN),
            size0 = box[0].hi - box[0].lo + 1,
            size1 = box[1].hi - box[1].lo + 1,
            minSize = MIN(size0, size1),
            maxSize = MAX(size0, size1)
          DO
            IF maxSize <= 2 OR minSize <= 1 THEN 
              FOR i0 := box[0].lo TO box[0].hi DO
                FOR i1 := box[1].lo TO box[1].hi DO
                  IF NOT InBox(i0, i1, orgBB)
                  AND NOT InBox(i0, i1, dstBB)
                  AND NOT seen[i0, i1] 
                  THEN
                    seen[i0, i1] := TRUE;
                    action(a, i0, i1, who[i0, i1])
                  END
                END
              END
            ELSE
              WITH
                p00m = Interpolate(p000, p001, 0.5),
                p0m1 = Interpolate(p001, p011, 0.5),
                pm11 = Interpolate(p011, p111, 0.5),
                p0mm = Interpolate(p00m, p0m1, 0.5),
                pmm1 = Interpolate(p0m1, pm11, 0.5),
                pmmm = Interpolate(p0mm, pmm1, 0.5)
              DO
                DoEnum(p000, p00m, p0mm, pmmm);
                DoEnum(pmmm, pmm1, pm11, p111);
              END
            END;
          END;
        END DoEnum;

      BEGIN
        DoEnum(p000, p001, p011, p111)
      END
    END
  END EnumBezierCells;
  
PROCEDURE InBox(i0, i1: INTEGER; READONLY box: IntBox): BOOL =
  BEGIN
    RETURN 
      i0 >= box[0].lo AND i0 <= box[0].hi AND 
      i1 >= box[1].lo AND i1 <= box[1].hi
  END InBox;
  
PROCEDURE Interpolate(READONLY p0, p1: RealPoint; t: REAL): RealPoint =
  BEGIN
    RETURN RealPoint{
      (1.0 - t)*p0[0] + t*p1[0],
      (1.0 - t)*p0[1] + t*p1[1]
    }
  END Interpolate;
  
PROCEDURE SegmentGridBox(READONLY p, q: RealPoint; READONLY gridN: Counts): IntBox =
  (* 
    Returns the ranges of *indices* of cells that contain the given segment. *)
  BEGIN
    WITH
      lo0 = MAX(0, FLOOR(MIN(p[0], q[0]))),
      lo1 = MAX(0, FLOOR(MIN(p[1], q[1]))),

      hi0 = MIN(gridN[0]-1, FLOOR(MAX(p[0], q[0]))),
      hi1 = MIN(gridN[1]-1, FLOOR(MAX(p[1], q[1])))
    DO
      RETURN IntBox{IntRange{lo0, hi0}, IntRange{lo1, hi1}}
    END
  END SegmentGridBox;

PROCEDURE BezierGridBox(READONLY p, q, r, s: RealPoint; READONLY gridN: Counts): IntBox =
  (*
    Returns the ranges of *indices* of cells that contain the given Bezier arc. *)
  BEGIN
    WITH
      lo0 = MAX(0, FLOOR(MIN(MIN(p[0], q[0]), MIN(r[0], s[0])))),
      lo1 = MAX(0, FLOOR(MIN(MIN(p[1], q[1]), MIN(r[1], s[1])))),

      hi0 = MIN(gridN[0]-1, FLOOR(MAX(MAX(p[0], q[0]), MAX(r[0], s[0])))),
      hi1 = MIN(gridN[1]-1, FLOOR(MAX(MAX(p[1], q[1]), MAX(r[1], s[1]))))
    DO
      RETURN IntBox{IntRange{lo0, hi0}, IntRange{lo1, hi1}}
    END
  END BezierGridBox;
  
PROCEDURE NodeArea(n: Node): NAT =
  BEGIN
    WITH
      ri = n.s[0].rad,
      ro = n.s[1].rad,
      ht = MAX(ri, ro)
    DO
      RETURN (ri + ro) * (ht + ht)
    END
  END NodeArea;

PROCEDURE Degree(n: Node; side: [0..1]): NAT =
  VAR deg: NAT := 0;
      a: Arc := n.s[side].arcs;
  BEGIN
    WHILE a # NIL DO 
      INC(deg);
      a := a.s[1 - side].next
    END;
    RETURN deg
  END Degree;

PROCEDURE MoveNode(<*UNUSED*> t: T; n: Node; gp: IntPoint) =
  BEGIN
    <* ASSERT NOT n.placed *>
    n.gp := gp
  END MoveNode;

PROCEDURE UnplaceNode(t: T; n: Node) =

  PROCEDURE RemoveNodeFromCell (* : NodeCellAction *) (
      n: Node;
      <*UNUSED*> i0: NAT; 
      <*UNUSED*> i1: NAT; 
      VAR cell: Fragment;
    ) =
    BEGIN
      <* ASSERT cell # NIL *>
      WITH 
        nc = cell.SubElement(n)
      DO
        <* ASSERT nc # cell *>
        cell := nc
      END
    END RemoveNodeFromCell;

  PROCEDURE RemoveArcFromCell (* : ArcCellAction *) (
      a: Arc;
      <*UNUSED*> i0: NAT;
      <*UNUSED*> i1: NAT; 
      VAR cell: Fragment;
    ) =
    BEGIN
      IF a = NIL THEN
        (* OK *)
      ELSE
        WITH nc = cell.SubElement(a) DO cell := nc END
      END
    END RemoveArcFromCell;

  <* FATAL Abort *>
  BEGIN
    <* ASSERT n.placed *>
    EnumNodeAndArcCells(t, n, RemoveNodeFromCell, RemoveArcFromCell, placedOnly := TRUE);
    IF NOT n.happy THEN RemoveNodeFromQueue(t.unhappy, n) END;
    n.happy := FALSE;
    n.placed := FALSE;
    PutNodeInQueue(t.unplaced, n);
    IF verbose OR NOT quiet THEN ReportUnplacement(n.label, n.gp) END
  END UnplaceNode;
  
PROCEDURE PlaceNode(t: T; n: Node; gp: IntPoint; happy: BOOLEAN; now: NAT) =

  PROCEDURE AddNodeToCell (* : NodeCellAction *) (
      n1: Node;
      <*UNUSED*> i0: NAT;
      <*UNUSED*> i1: NAT; 
      VAR cell: Fragment;
    ) =
    BEGIN
      <* ASSERT n1 = n *> (* To force AddArcToCell to be a closure *)
      IF cell = NIL THEN
        cell := n
      ELSE
        WITH 
          nc = cell.AddElement(n)
        DO
          <* ASSERT nc # cell *>
          cell := nc
        END
      END
    END AddNodeToCell;

  PROCEDURE AddArcToCell (* : ArcCellAction *) (
      a: Arc;
      <*UNUSED*> i0: NAT;
      <*UNUSED*> i1: NAT; 
      VAR cell: Fragment;
    ) =
    BEGIN
      <* ASSERT a.s[0].node = n OR a.s[1].node = n *> (* To force AddArcToCell to be a closure *)
      IF cell = NIL THEN
        cell := a
      ELSE
        WITH nc = cell.AddElement(a) DO cell := nc END
      END
    END AddArcToCell;

  <* FATAL Abort *>
  BEGIN
    <* ASSERT NOT n.placed *>
    <* ASSERT NOT n.happy *>
    <* ASSERT gp[0] >= n.range[0].lo AND gp[0] <= n.range[0].hi *>
    <* ASSERT gp[1] >= n.range[1].lo AND gp[1] <= n.range[1].hi *>
    n.gp := gp;
    UnplaceConflictingNodes(t, n);
    n.placed := TRUE;
    n.happy := happy;
    n.when := now;
    EnumNodeAndArcCells(t, n, AddNodeToCell, AddArcToCell, placedOnly := TRUE);
    RemoveNodeFromQueue(t.unplaced, n);
    IF NOT n.happy THEN PutNodeInQueue(t.unhappy, n) END;
    IF verbose OR NOT quiet THEN ReportPlacement(n.label, n.gp) END
  END PlaceNode;
  
PROCEDURE UnplaceConflictingNodes(t: T; n: Node) =
  (*
    Unplaces all placed nodes that collide with the currently
    unplaced node "n", and puts them in the unplaced node queue: *)

  PROCEDURE ProcessOtherNode(m: Node) =
    BEGIN
      <* ASSERT m # n *>
      <* ASSERT m.placed *>
      t.UnplaceNode(m)
    END ProcessOtherNode;
  
  PROCEDURE ProcessCell (* : NodeCellAction *) (
      n1: Node;
      <*UNUSED*> i0: NAT;
      <*UNUSED*> i1: NAT; 
      VAR cell: Fragment;
    ) =
    BEGIN
      <* ASSERT n = n1 *>
      IF cell # NIL THEN
        cell.EnumElements(arcAction := NIL, nodeAction := ProcessOtherNode)
      END
    END ProcessCell;

  <* FATAL Abort *>
  BEGIN
    <* ASSERT NOT n.placed *>
    t.EnumNodeCells(n, action := ProcessCell);
  END UnplaceConflictingNodes;

PROCEDURE EnumNodeAndArcCells(
    t: T; 
    n: Node; 
    nodeAction: NodeCellAction;
    arcAction: ArcCellAction;
    placedOnly: BOOL := FALSE;
  ) RAISES {Abort} =
  (*
    Enumerates all cells covered by the given node, and
    by the arcs incident to it.  If "placedOnly" is true, 
    considers only arcs  whose endpoints are both placed. *)
  BEGIN
    IF nodeAction # NIL THEN EnumNodeCells(t, n, nodeAction) END;
    IF arcAction # NIL THEN
      FOR side := 0 TO 1 DO 
        VAR a: Arc := n.s[side].arcs;
        BEGIN
          WHILE a # NIL DO
            IF (NOT placedOnly) OR ArcIsPlaced(a) THEN 
              EnumArcCells(t, a, arcAction)
            END;
            a := a.s[1-side].next
          END
        END;
      END;
    END
  END EnumNodeAndArcCells;
  
CONST DebugEnergy = TRUE;

PROCEDURE ComputeIdealPosition(t: T; n: Node; READONLY range: IntBox): IntPoint =
  VAR 
    energy, xForce, yForce: REAL;
    energyBest, xForceBest, yForceBest: REAL;
    gp, gpBest: IntPoint;
    trials: NAT := 0;
    estMinEnergy: REAL := 0.0;
    xStep, yStep, dMin, dMax: REAL;
  CONST MaxTrials = 5;
  BEGIN
    <* ASSERT NOT n.placed *>
    energyBest := LAST(REAL);
    gpBest := Nowhere;
    gp := n.gp;
    dMin := 1.0/FLOAT(MAX(t.gridN[0], t.gridN[1]));
    dMax := dMin;
    REPEAT
      (* Compute energy and force at trial position "gp": *)
      ComputeEnergyAndForce(n, gp, t.gridN, TRUE, energy, xForce, yForce, cutoff := LAST(REAL));
      INC(trials);
      IF energy < energyBest THEN 
        IF gpBest = Nowhere THEN
          estMinEnergy := 0.0
        ELSE
          estMinEnergy := energy - (energyBest - energy);
        END;
        WITH
          f2 = xForce*xForce + yForce*yForce,
          s = (energy - estMinEnergy)/f2
        DO
          xStep := 2.0 * xForce * s;
          yStep := 2.0 * yForce * s;
          (* Clip step length to dMax: *)
          WITH
            step2 = xStep*xStep + yStep*yStep
          DO
            IF step2 > dMax * dMax THEN 
              WITH r = dMax/Sqrt(step2) DO 
                xStep := xStep * r;
                yStep := yStep * r
              END;
              dMax := 2.0 * dMax
            ELSE
              dMax := MAX(dMin, 2.0 * Sqrt(step2))
            END
          END
        END;
        energyBest := energy;
        xForceBest := xForce;
        yForceBest := yForce;
        gpBest := gp
      ELSE
        xStep := 0.5 * xStep;
        yStep := 0.5 * yStep;
        dMax := MAX(dMin, 0.5 * dMax)
      END;
      gp[0] := gpBest[0] + ROUND(xStep * FLOAT(t.gridN[0]));
      gp[1] := gpBest[1] + ROUND(yStep * FLOAT(t.gridN[1]));
      
      gp := ClipIntPoint(gp, range);
    UNTIL gp = gpBest OR trials >= MaxTrials;
    RETURN gpBest
  END ComputeIdealPosition;
  
CONST 
  EnableCutoff = FALSE; (* Enable energy "cutoff" optimizations *)

PROCEDURE ComputeEnergyAndForce(
    n: Node;
    gp: IntPoint;
    READONLY gridN: Counts;
    dirFixed: BOOL;
    VAR (*OUT*) energy: REAL; 
    VAR (*OUT*) xForce: REAL;
    VAR (*OUT*) yForce: REAL;
    cutoff: REAL;
  ) =
  (*
    Computes the total potential energy for node "n" if it were plaed at "gp".
    The energy includes the node's gravitational energy and half
    of the strain energy of the arcs incident to "n". 
    
    If "dirFixed", assumes that arc directions at "n" are fixed;
    else assumes the arc curvature there is zero.
    
    May abort the computation, with a partial energy amount, 
    if it finds out that the total energy would be bigger than "cutoff".
    *)

  PROCEDURE DoComputeEnergy() =
    VAR a: Arc;
    BEGIN
      energy := 0.0;
      xForce := 0.0;
      yForce := 0.0;
      ComputeGravity(gp, gridN, energy, xForce, yForce);
      IF EnableCutoff AND energy >= cutoff THEN RETURN END;
      FOR side := 0 TO 1 DO
        a := n.s[side].arcs;
        WHILE a # NIL DO
          <* ASSERT a.s[1-side].node = n *>
          WITH
            m = a.s[side].node,
            mPos = m.gp,
            mDir = a.s[side].dir + m.tilt,
            nPos = gp, 
            nDir = a.s[1-side].dir + n.tilt
          DO
            ComputeArcStrain(
              a,
              mPos, mDir, nPos, nDir, 
              gridN, dirFixed, 
              energy, xForce, yForce
            );
          END; 
          IF EnableCutoff AND energy >= cutoff THEN RETURN END;
          a := a.s[1-side].next
        END;
      END;
    END DoComputeEnergy;

  BEGIN
    DoComputeEnergy();
    IF DebugEnergy AND verbose THEN
      ReportEnergyComputation(n.label, gp, dirFixed, energy, xForce, yForce)
    END;
  END ComputeEnergyAndForce;

PROCEDURE ComputeGravity(
    READONLY gp: IntPoint; 
    READONLY gridN: Counts;
    VAR (*INOUT*) energy: REAL; 
    VAR (*INOUT*) xForce: REAL;
    VAR (*INOUT*) yForce: REAL;
  ) =
  (* 
    Computes potential energy and force for placing a node at "gp",
    ignoring all other nodes and arcs.  Adds results to "energy",
    "xForce", "yForce". *)
  CONST GravityModulus = 0.01;  (* Force constant for displacement from center of grid. *)
  BEGIN
    WITH
      dx = FLOAT(gp[0])/FLOAT(gridN[0]) - 0.5,
      dy = FLOAT(gp[1])/FLOAT(gridN[1]) - 0.5
    DO
      energy := energy + GravityModulus * (dx * dx + dy * dy);
      xForce := xForce - 2.0 * GravityModulus * dx;
      yForce := yForce - 2.0 * GravityModulus * dy;
    END;
  END ComputeGravity;  

PROCEDURE ComputeArcStrain(
    a: Arc;
    READONLY mPos: IntPoint; 
    mDir: REAL;
    READONLY nPos: IntPoint; 
    nDir: REAL;
    READONLY gridN: Counts;
    dirFixed: BOOL;
    VAR (*INOUT*) energy: REAL;
    VAR (*INOUT*) xForce: REAL;
    VAR (*INOUT*) yForce: REAL;
  ) =
  (*
    Computes the strain energy of an arc that connects a fixed node
    "m" at "mPos" with a movable node "n" at "nPos".  
    
    Computes also the derivatives of the strain with repect to the
    position of the moveable node "n".  (Coordinates are scaled so
    that [0 __ 1] spans the whole window.)
    
    Adds results to "energy", "xForce", "yForce".  
    
    Assumes that the direction of the arc at node "m" is fixed at "mDir".
    If "dirFixed", assumes the direction at node "n" is "nDir"; else
    ignores "nDir", and computes the direction at "n" assuming that the 
    curvature there is zero. *)

  CONST StretchModulus = 10.0;
        TiltModulus = 0.1;
        BendModulus = 0.05;
  VAR sense: REAL;
      angn, dangn_dxn, dangn_dyn: REAL;
  BEGIN
    IF mPos[0] < nPos[0] THEN sense := +1.0 ELSE sense := -1.0 END;
    WITH
      xm = FLOAT(mPos[0])/FLOAT(gridN[0]),
      ym = FLOAT(mPos[1])/FLOAT(gridN[1]),
      
      xn = FLOAT(nPos[0])/FLOAT(gridN[0]),
      yn = FLOAT(nPos[1])/FLOAT(gridN[1]),
      
      (* Vector from "m" to "n": *)
      rx = xn - xm,
      ry = yn - ym,
      
      eps2 = Sqr(1.0/FLOAT(gridN[0])) + Sqr(1.0/FLOAT(gridN[1])),
      
      r2 = rx * rx + ry * ry + eps2,    
      dr2_dxn = 2.0 * rx,       
      dr2_dyn = 2.0 * ry,
      
      r = Sqrt(r2),
      dr_dxn = rx/r,
      dr_dyn = ry/r,
      
      (* "angs" is the direction of the straight line vector
        from the *left* endpoint to the *right* endpoint: *)
      angs = Dir(sense * rx, sense * ry),
      dangs_dxn = - sense * Sin(angs)/r,
      dangs_dyn = sense * Cos(angs)/r,

      (* "angm" is the local arc direction at the fixed endpoint: *)
      angm = mDir,
      dangm_dxn = 0.0,
      dangm_dyn = 0.0

    DO
      
      (* "angn" is the local arc direction at the variable endpoint: *)
      IF dirFixed THEN
        angn := nDir;
        dangn_dxn := 0.0;
        dangn_dyn := 0.0
      ELSE
        angn := angs - 0.5*(angm - angs);
        dangn_dxn := dangs_dxn - 0.5*(dangm_dxn - dangs_dxn);
        dangn_dyn := dangs_dyn - 0.5*(dangm_dyn - dangs_dyn);
      END;
      
      WITH
        (* "SE" is the stretching energy: *)
        SE = StretchModulus * r2,
        dSE_dxn = StretchModulus * dr2_dxn,
        dSE_dyn = StretchModulus * dr2_dyn,
        
        (* "TE" is the tilt energy: *)
        TE = TiltModulus * angs * angs,
        dTE_dxn = 2.0 * TiltModulus * angs * dangs_dxn,
        dTE_dyn = 2.0 * TiltModulus * angs * dangs_dyn,

        (* "turn" is the total slope variation: *)
        turn = (angn - angm)/2.0,
        dturn_dxn = (dangn_dxn - dangm_dxn)/2.0,
        dturn_dyn = (dangn_dyn - dangm_dyn)/2.0,

        (* "wave" is the total wawiness: *)
        wave = (angn + angm)/2.0 - angs,
        dwave_dxn = (dangn_dxn + dangm_dxn)/2.0 - dangs_dxn,
        dwave_dyn = (dangn_dyn + dangm_dyn)/2.0 - dangs_dyn,

        (* "w" is the total squared bending: *)
        w = Sqr(turn) + 3.0 * Sqr(wave),
        dw_dxn = 2.0 * turn * dturn_dxn + 6.0 * wave * dwave_dxn,
        dw_dyn = 2.0 * turn * dturn_dyn + 6.0 * wave * dwave_dyn,
        
        (* "BE" is the bending energy: *)
        BE = 4.0 * BendModulus * w / r,
        dBE_dxn = 4.0 * BendModulus * (dw_dxn - w * dr_dxn)/r,
        dBE_dyn = 4.0 * BendModulus * (dw_dyn - w * dr_dyn)/r

      DO
        IF verbose AND DebugEnergy THEN 
          ReportArcStrain(a, 
            rx, ry, angm, angn,
            SE, dSE_dxn, dSE_dyn, 
            TE, dTE_dxn, dTE_dyn,
            BE, dBE_dxn, dBE_dyn
          )
        END;
        energy := energy + (SE + TE + BE);
        xForce := xForce - (dSE_dxn + dTE_dxn + dBE_dxn);
        yForce := yForce - (dSE_dyn + dTE_dyn + dBE_dyn);
      END;
    END    
  END ComputeArcStrain;

PROCEDURE PositionPenalty(t: T; n: Node; cutoff: REAL): REAL =
  VAR energy, xForce, yForce: REAL;
  BEGIN
    ComputeEnergyAndForce(n, n.gp, t.gridN, TRUE, energy, xForce, yForce, cutoff);
    RETURN energy
  END PositionPenalty;

CONST
  NodeNodeCollisionPenalty = 2000.0; (* Basic penalty for node/node overlap *)
  ArcNodeCollisionPenalty = 1000.0;  (* Penalty for node/arc overlap *)
  ArcArcCollisionPenalty = 100.0;    (* Penalty for arc/arc overlap *)

PROCEDURE CollisionPenalty(t: T; n: Node; now: NAT; cutoff: REAL): REAL =
  (*
    Penalty (due to colisions with already placed nodes)
    if we were to place node /n/ at position /gp/. *)
  VAR
    nnHits: ARRAY [0..10] OF Node;
    nnHitCount: NAT := 0;
    naHits: ARRAY [0..15] OF Arc;
    naHitCount: NAT := 0;
    anHits: ARRAY [0..15] OF Node;
    anHitCount: NAT := 0;
    aaHits: ARRAY [0..30] OF Arc;
    aaHitCells: ARRAY [0..30] OF NAT;
    aaHitCount: NAT := 0;
    penalty: REAL := 0.0;
    
  PROCEDURE ProcessNodeCell (
      n1: Node; 
      <*UNUSED*> i0: NAT;
      <*UNUSED*> i1: NAT; 
      VAR cell: Fragment;
    ) RAISES {Abort} =
    
    PROCEDURE ProcessNodeOverNode(m: Node) RAISES {Abort} =
      BEGIN
        <* ASSERT m # n *>
        FOR i := 0 TO nnHitCount-1 DO 
          IF nnHits[i] = m THEN RETURN END;
        END;
        IF nnHitCount = NUMBER(nnHits) THEN RAISE Abort END;
        nnHits[nnHitCount] := m;
        INC(nnHitCount);
        penalty := penalty + ComputeNodeNodeCollisionPenalty(t, n, m, now);
        IF EnableCutoff AND penalty > cutoff THEN RAISE Abort END
      END ProcessNodeOverNode;
      
    PROCEDURE ProcessNodeOverArc(b: Arc) RAISES {Abort} =
      BEGIN
        <* ASSERT b.s[0].node # n AND b.s[1].node # n *>
        FOR i := 0 TO naHitCount-1 DO 
          IF naHits[i] = b THEN RETURN END;
        END;
        IF naHitCount = NUMBER(naHits) THEN RAISE Abort END;
        naHits[naHitCount] := b;
        INC(naHitCount);
        penalty := penalty + ArcNodeCollisionPenalty;
        IF EnableCutoff AND penalty > cutoff THEN RAISE Abort END
      END ProcessNodeOverArc;

    BEGIN
      <* ASSERT n1 = n *>
      IF cell # NIL THEN
        cell.EnumElements(nodeAction := ProcessNodeOverNode, arcAction := ProcessNodeOverArc)
      END;
    END ProcessNodeCell;
    
  PROCEDURE ProcessArcCell (
      a: Arc; 
      (*UNUSED*) i0: NAT;
      (*UNUSED*) i1: NAT;
      VAR cell: Fragment;
    ) RAISES {Abort} =

    PROCEDURE ProcessArcOverNode(m: Node) RAISES {Abort} =
      BEGIN
        <* ASSERT a.s[0].node # m AND a.s[1].node # m *>
        FOR i := 0 TO anHitCount-1 DO 
          IF anHits[i] = m THEN RETURN END;
        END;
        IF anHitCount = NUMBER(anHits) THEN RAISE Abort END;
        anHits[anHitCount] := m;
        INC(anHitCount);
        penalty := penalty + ArcNodeCollisionPenalty;
        IF EnableCutoff AND penalty > cutoff THEN RAISE Abort END
      END ProcessArcOverNode;
      
    PROCEDURE ProcessArcOverArc(b: Arc) RAISES {Abort} =
      BEGIN
        <* ASSERT b # NIL *>
        <* ASSERT a # b *>
        IF verbose THEN ReportArcArcCollision(a, b, i0, i1) END;
        FOR i := 0 TO aaHitCount-1 DO 
          IF aaHits[i] = b THEN 
            INC(aaHitCells[i]);
            IF aaHitCells[i] > 4 THEN
              penalty := penalty + ArcArcCollisionPenalty;
              IF EnableCutoff AND penalty > cutoff THEN RAISE Abort END;
            END;
            RETURN
          END;
        END;
        IF aaHitCount = NUMBER(aaHits) THEN RAISE Abort END;
        aaHits[aaHitCount] := b;
        INC(aaHitCount);
        penalty := penalty + ArcArcCollisionPenalty;
        IF EnableCutoff AND penalty > cutoff THEN RAISE Abort END
      END ProcessArcOverArc;

    BEGIN
      IF cell # NIL THEN
        cell.EnumElements(arcAction := ProcessArcOverArc, nodeAction := ProcessArcOverNode)
      END;
    END ProcessArcCell;
    
  BEGIN
    <* ASSERT NOT n.placed *>
    TRY
      t.EnumNodeCells(n, action := ProcessNodeCell);
      VAR a: Arc;
      BEGIN
        FOR side := 0 TO 1 DO 
          a := n.s[side].arcs;
          WHILE a # NIL DO
            <* ASSERT a.s[1-side].node = n *>
            IF a.s[side].node.placed THEN
              aaHitCount := 0;
              anHitCount := 0;
              FOR i := 0 TO LAST(aaHitCells) DO aaHitCells[i] := 0 END;
              t.EnumArcCells(a, action := ProcessArcCell)
            END;
            a := a.s[1-side].next
          END
        END
      END;
      RETURN penalty
    EXCEPT
      Abort => RETURN ExcessivePenalty
    END;
  END CollisionPenalty;

PROCEDURE ComputeNodeNodeCollisionPenalty(t: T; n, m: Node; now: NAT): REAL =
  (*
    Penalty for placing a node "n" where it would cover one or more cells 
    of another node "m" (thus forcing the latter to be unplaced).
    Returns ExcessivePenalty if "m" is permanently fixed.
    (Should do the same if "m" can't be displaced far enough...)
    *)
  BEGIN
    <* ASSERT m # NIL *>
    <* ASSERT n # NIL *>
    <* ASSERT m.placed *>
    <* ASSERT NOT n.placed *>
    WITH
      nbb = NodeBBox(t, n),
      mri = m.s[0].rad,
      mro = m.s[1].rad,
      mht = m.ht
    DO
      IF  m.range[0].lo + mro - 1 >= nbb[0].lo 
      AND m.range[0].hi - mri <= nbb[0].hi
      AND m.range[1].lo + mht - 1 >= nbb[1].lo
      AND m.range[1].hi - mht <= nbb[1].hi
      THEN
        (* Node "m" has nowhere else to go, no way: *)
        RETURN ExcessivePenalty
      ELSE
        <* ASSERT m.range[0].lo < m.range[0].hi OR m.range[1].lo < m.range[1].hi *>
        RETURN NodeNodeCollisionPenalty*FLOAT(m.when)/FLOAT(MAX(1, now))
      END
    END
  END ComputeNodeNodeCollisionPenalty;
  
PROCEDURE ComputeArcControlPoints(a: Arc): ArcControlPoints =
  BEGIN
    WITH
      org = a.s[0].node,
      dst = a.s[1].node,
      
      dir0 = a.s[0].dir + org.tilt,
      dir1 = a.s[1].dir + dst.tilt,
      
      r0 = FLOAT(org.s[1].rad),
      r1 = FLOAT(dst.s[0].rad),
      
      x0 = FLOAT(org.gp[0]),
      y0 = FLOAT(org.gp[1]),
      
      x1 = FLOAT(dst.gp[0]),
      y1 = FLOAT(dst.gp[1]),

      xa = x0 + r0 * Cos(dir0),
      ya = y0 + r0 * Sin(dir0),

      xd = x1,
      yd = y1,

      xb = (2.0*xa + xd)/3.0,
      yb = ya + (xb - xa)*Tan(dir0),

      xc = (xa + 2.0*xd)/3.0,
      yc = yd - (xd - xc)*Tan(dir1)

    DO
      RETURN ArcControlPoints{
        RealPoint{x0, y0},
        RealPoint{xa, ya},
        RealPoint{xb, yb},
        RealPoint{xc, yc},
        RealPoint{xd, yd}
      }
    END;
  END ComputeArcControlPoints;

PROCEDURE ArcIsPlaced(a: Arc): BOOL =
  BEGIN
    RETURN a.s[0].node.placed AND a.s[1].node.placed
  END ArcIsPlaced;
  
PROCEDURE MakeNodeUnhappy(t: T; n: Node) =
  BEGIN
    <* ASSERT n.placed *>
    IF n.happy THEN
      n.happy := FALSE;
      PutNodeInQueue(t.unhappy, n);
      IF verbose THEN ReportUnhappiness(n.label) END;
    END
  END MakeNodeUnhappy;

(* NODE QUEUES *)

PROCEDURE NewQueue(maxNodes: NAT): NodeQueue =
  BEGIN
    WITH queue = NEW(NodeQueue,
        node := NEW(REF ARRAY OF Node, maxNodes + 2),
        front := 0,
        rear := 0
      )
    DO 
      RETURN queue
    END
  END NewQueue;

PROCEDURE QueueLength(queue: NodeQueue): NAT =
  BEGIN
    WITH
      size = NUMBER(queue.node^),
      n = (queue.rear + size - queue.front) MOD size
    DO
      RETURN n
    END
  END QueueLength;

PROCEDURE FirstNodeInQueue(queue: NodeQueue): Node =
  BEGIN
    IF queue.front = queue.rear THEN
      RETURN NIL
    ELSE
      RETURN queue.node[queue.front]
    END
  END FirstNodeInQueue;

PROCEDURE RemoveNodeFromQueue(queue: NodeQueue; n: Node) =
  BEGIN
    <* ASSERT n # NIL *>
    <* ASSERT queue.front # queue.rear *>
    WITH
      node = queue.node^,
      size = NUMBER(node)
    DO
      VAR j: NAT := queue.front;
          m: Node := node[j];
          t: Node;
      BEGIN
        WHILE m # n DO
          INC(j);
          IF j >= size THEN j := 0 END;
          <* ASSERT j # queue.rear *>
          t := node[j]; node[j] := m; m := t;
        END;
        INC(queue.front);
        IF queue.front >= size THEN queue.front := 0 END;
      END
    END
  END RemoveNodeFromQueue;

PROCEDURE PutNodeInQueue(queue: NodeQueue; n: Node) =
  BEGIN
    <* ASSERT n # NIL *>
    queue.node[queue.rear] := n;
    INC(queue.rear);
    IF queue.rear > LAST(queue.node^) THEN queue.rear := 0 END;
    <* ASSERT queue.rear # queue.front *>
  END PutNodeInQueue;

PROCEDURE CycleNodesInQueue(queue: NodeQueue) =
  VAR n: Node;
  BEGIN
    IF queue.rear # queue.front THEN
      n := queue.node[queue.front];
      <* ASSERT n # NIL *>
      INC(queue.front);
      IF queue.front > LAST(queue.node^) THEN queue.front := 0 END;
      queue.node[queue.rear] := n;
      INC(queue.rear);
      IF queue.rear > LAST(queue.node^) THEN queue.rear := 0 END;
      <* ASSERT queue.rear # queue.front *>
    END;
  END CycleNodesInQueue;
  
BEGIN
END ADDiagram.

