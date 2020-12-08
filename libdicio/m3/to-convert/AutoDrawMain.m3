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

(* Last modified on Thu Aug 20 17:01:18 PDT 1992 by stolfi                  *)

MODULE AutoDrawMain EXPORTS Main;

IMPORT
  ADReport, ADDiagram,
  ADPlot,
  Basics, Reduced, Encoding, PlainEncoding, 
  Rd, Wr, Text, Fmt, FileStream, ParseParams, 
  Random, Thread, Scan;
FROM Basics IMPORT INT, NAT, BOOL;
FROM Reduced IMPORT State, NullState;
FROM Stdio IMPORT stdout, stderr;
FROM ADTypes IMPORT
  Dims, Counts, IntPoint, IntRange, IntBox, Nowhere, 
  RealPoint, ArcControlPoints, ClipIntPoint;
FROM ADDiagram IMPORT
  Node, Arc, QueueLength,
  ComputeIdealPosition, PositionPenalty, CollisionPenalty;
FROM ADMath IMPORT 
  Cos, Sin, Dir, L;
FROM ADReport IMPORT
  ReportQueueStatus,
  ReportSuccess, ReportAttemptToPlace, ReportTentativePlacement,
  ReportIdealPosition, ReportCurrentRange, 
  ReportInsufficientHorRange, ReportFailureToPlace, 
  ErrorNotEnoughSpace, ErrorNotEnoughGridCells,
  PrintGridData;

TYPE
  Options = RECORD
      autName: TEXT;        (* Name of automaton *)
      verbose: BOOL;        (* TRUE to mumble a lot while thinking. *)
      quiet: BOOL;          (* TRUE to be very quiet while thinking. *)
      gridDims: Dims;       (* Size of drawing area, in meters *)
      writeCoords: BOOL;    (* TRUE to write node/arc coordinate files  *)
      gridStep: REAL;       (* Mesh of grid of allowed state positions *)
      nodeRadius: REAL;     (* Radius of node *)
      arcSpacing: REAL;     (* Spacing between arcs at the output radius *)
      iterations: NAT;      (* Nominal max number of iterations *)
      showStart: NAT;       (* Number of iterations to first snapshot *)
      showStep: NAT;        (* Number of iterations between snapshots *)
      showCount: NAT;       (* Number of snapshots (excluding final one) *)
    END;

VAR o: Options;

CONST
  MaxFanWidth = 0.90 * 3.1415926;  (* Max angular spread of input/output edges, in radians *)
  Degree = 3.1415926/180.0;        (* One degree, in radians *)

PROCEDURE Main() =
  <* FATAL Wr.Failure, Rd.Failure, Thread.Alerted *>
  BEGIN
    o := GetOptions();
    ADDiagram.verbose := o.verbose;
    ADDiagram.quiet := o.quiet;
    WITH
      t = BuildDiagram(),
      psFile = OpenPSFile()
    DO
      ComputeNodePositionsAndDraw(t, psFile);
      IF o.writeCoords THEN WriteCoordFiles(t) END;
    END;
    Wr.Flush(stderr);
    Wr.Flush(stdout);
  END Main;
  
PROCEDURE BuildDiagram(): ADDiagram.T =
  <* FATAL Wr.Failure, Rd.Failure, Thread.Alerted *>
  VAR aut: Reduced.T;
  BEGIN

    (* Load the automaton: *)
    <* ASSERT NOT Text.Empty(o.autName) *>
    WITH
      loadFileName = o.autName & ".dmp"
    DO
      Wr.PutText(stderr, "\n=== loading automaton from " & loadFileName & " ===\n\n");
      WITH rd = FileStream.OpenRead(loadFileName) DO
        aut := Reduced.Load(rd)
      END;
    END;

    WITH 
      root = aut.Root(),
      ct = aut.Count(ARRAY OF State{root}) 
    DO

      (* Print report *)
      Wr.PutText(stderr, "root state = " & Fmt.Pad(Fmt.Int(aut.Root()), 8) & "\n");
      Wr.PutText(stderr, "states =     " & Fmt.Pad(Fmt.Int(ct.states), 8) & "\n");
      Wr.PutText(stderr, "finals =     " & Fmt.Pad(Fmt.Int(ct.finals), 8) & "\n");
      Wr.PutText(stderr, "arcs =       " & Fmt.Pad(Fmt.Int(ct.arcs), 8) & "\n");

      (* Convert automaton into diagram: *)
      WITH
        gridN = Counts{
          ROUND(o.gridDims[0]/o.gridStep),
          ROUND(o.gridDims[1]/o.gridStep)
        },
        t = ADDiagram.New(maxLabel := root, maxNodes := ct.states, gridN := gridN)
      DO
        CollectReachableStates(t, aut);
        RETURN t
      END
    END;
  END BuildDiagram;

PROCEDURE CollectReachableStates(t: ADDiagram.T; aut: Reduced.T) =
  (*
    Enumerates all true states and arcs of the automaton /aut/
    reachable from /aut.Root()/, and their arcs, and adds them to /t/. *)

  PROCEDURE InRadius(<*UNUSED*> indeg: NAT): REAL =
    BEGIN
      RETURN 1.5 * o.nodeRadius
    END InRadius;

  PROCEDURE OutRadius(outdeg: NAT): REAL =
    BEGIN
      RETURN MAX(
        o.nodeRadius + o.arcSpacing, 
        (FLOAT(outdeg + 1) * o.arcSpacing) / MaxFanWidth
      )
    END OutRadius;
    
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    IF o.verbose THEN Wr.PutText(stderr, "Collecting node and arc data...\n") END;
    WITH
      r = aut.Root(),
      encoding = PlainEncoding.New()
    DO
      FOR s := 1 TO r DO
        IF aut.NPrefs(s) > 0 THEN
          WITH
            ri = CEILING(InRadius(aut.InDeg(s))/o.gridStep),
            ro = CEILING(OutRadius(aut.OutDeg(s))/o.gridStep),
            ht = MAX(ri, ro),
            n = t.AddNode(s, aut.Final(s), ri, ro, ht)
          DO
            AddOutputArcs(t, aut, encoding, s)
          END
        END
      END
    END
  END CollectReachableStates;

PROCEDURE AddOutputArcs(
    t: ADDiagram.T; 
    aut: Reduced.T; 
    encoding: PlainEncoding.T; 
    s: Reduced.State;
  ) =

  PROCEDURE ProcessSucc (* : Reduced.ArcAction *) (<*UNUSED*> i: NAT; a: Reduced.Arc) =
    <* FATAL Encoding.BadLetter *>
    BEGIN
      <* ASSERT a.dest # NullState *>
      WITH 
        label = encoding.LetterToChar(a.letter),
        org = t.node[s],
        dst = t.node[a.dest]
      DO
        <* ASSERT org # NIL *>
        <* ASSERT dst # NIL *>
        EVAL t.AddArc(label, org, dst)
      END
    END ProcessSucc;

  <* FATAL Basics.Skip, Basics.Abort *>
  BEGIN
    <* ASSERT s # NullState *>
    aut.EnumOutArcs(s, action := ProcessSucc);
  END AddOutputArcs;

PROCEDURE ComputeInitialPositions(t: ADDiagram.T; now: NAT) =
  (*
    Computes initial positions for all nodes; places and fixes the sources and sinks, 
    leaves rest unplaced. *)
  VAR nSources: NAT := 0;  (* Source (non-sink) nodes *)
      nSinks: NAT := 0;    (* Sink (non-source) nodes *)
      nSpouts: NAT := 0;   (* Sink-AND-source nodes *)
      wdSpouts: NAT := 0;  (* Total width of spouts *)
      htSources: NAT := 0; (* Total height of sources *)
      htSinks: NAT := 0;   (* Total height of sinks *)
  BEGIN
    WITH 
      node = t.node^,
      gridN = t.gridN
    DO 
      (* Counts sources, sinks, and spouts: *)
      FOR s := 1 TO LAST(node) DO 
        WITH n = node[s] DO 
          IF n # NIL THEN 
            IF n.s[0].arcs = NIL AND n.s[1].arcs = NIL THEN
              INC(nSpouts);
              INC(wdSpouts, n.s[0].rad + 1 + n.s[1].rad)
            ELSIF n.s[0].arcs = NIL THEN
              INC(nSources);
              INC(htSources, n.ht + 1 + n.ht)
            ELSIF n.s[1].arcs = NIL THEN 
              INC(nSinks);
              INC(htSinks, n.ht + 1 + n.ht)
            END
          END
        END
      END;
      
      (* Places sources, sinks, spouts: *)
      VAR gp: IntPoint; (* Initial position *)
          fix: ARRAY [0..1] OF BOOL;   (* TRUE = fix coordinate *)
          place: BOOL;  (* TRUE = place node, FALSE = just set gp *)
          
          dxSpout: INT := MAX(0, gridN[0] - wdSpouts) DIV (nSpouts + 1);
          dySource: INT := MAX(0, gridN[1] - htSources) DIV (nSources + 1);
          dySink: INT := MAX(0, gridN[1] - htSinks) DIV (nSinks + 1);
          
          xSpout: INT := dxSpout;
          ySource: INT := gridN[1] - dySource;
          ySink: INT := gridN[1] - dySink;
      BEGIN
        FOR s := 1 TO LAST(node) DO
          WITH n = node[s] DO
            IF n # NIL THEN
              <* ASSERT NOT n.placed *>
              (* Decide position /gp/ and whether it is fixed or mobile: *)
              IF n.s[0].arcs = NIL AND n.s[1].arcs = NIL THEN
                (* Spout node, place along bottom edge: *)
                gp := IntPoint{xSpout + n.s[0].rad, n.ht};
                xSpout := xSpout + n.s[0].rad + 1 + n.s[1].rad + dxSpout;
                place := TRUE;
                fix[0] := TRUE;
                fix[1] := TRUE;
              ELSIF n.s[0].arcs = NIL THEN
                (* Source node, place along left edge: *)
                gp := IntPoint{n.s[0].rad, ySource - n.ht - 1};
                ySource := ySource - n.ht - 1 - n.ht - dySource;
                place := TRUE;
                fix[0] := TRUE;
                fix[1] := (nSources = 1);
              ELSIF n.s[1].arcs = NIL THEN
                (* Sink node, place along right edge: *)
                gp := IntPoint{gridN[0] - n.s[1].rad, ySink - n.ht - 1};
                ySink := ySink - n.ht - 1 - n.ht - dySink;
                place := TRUE;
                fix[0] := TRUE;
                fix[1] := (nSinks = 1);
              ELSE
                gp := IntPoint{
                  (n.range[0].lo + n.range[0].hi) DIV 2,
                  (n.range[1].lo + n.range[1].hi) DIV 2
                };
                place := FALSE;
                fix[0] := FALSE;
                fix[1] := FALSE;
              END;
              
              (* Check if computed gp lies in valid range: *)
              FOR axis := 0 TO 1 DO
                IF gp[axis] < n.range[axis].lo OR gp[axis] > n.range[axis].hi THEN
                  gp[axis] := (n.range[axis].lo + n.range[axis].hi) DIV 2;
                  place := FALSE;
                  fix[0] := FALSE;
                  fix[1] := FALSE;
                END;
              END;

              (* Now place node, if appropriate: *)
              IF gp # Nowhere THEN
                IF place THEN
                  FOR axis := 0 TO 1 DO
                    IF fix[axis] THEN
                      n.range[axis] := IntRange{gp[axis], gp[axis]}
                    END
                  END;
                  t.PlaceNode(n, gp, happy := FALSE, now := now)
                ELSE
                  t.MoveNode(n, gp)
                END;
              END;
            END;
          END
        END
      END;
    END
  END ComputeInitialPositions;
  
PROCEDURE BasicArcDir(i, deg: NAT): REAL =
  (*
    Computes the "basic"  direction for the "i"th incoming/outgoing
    arc out of "n" such arcs.  Does not take into account node coordinates,
    and assumes the node's "tilt" is 0.0. *)
  BEGIN
    <* ASSERT 0 <= i AND i < deg *>
    WITH
      frac = FLOAT(i+1)/FLOAT(deg+1),
      angle = (frac - 0.5) * MaxFanWidth
    DO
      RETURN angle
    END
  END BasicArcDir;

PROCEDURE ComputeInitialArcDirs(t: ADDiagram.T) =

  PROCEDURE ProcessNode(n: Node; side: [0..1]) =
  
    VAR deg: NAT := ADDiagram.Degree(n, side);
    
    PROCEDURE ProcessArc (a: Arc; i: NAT) =
      (*
        Processes the /i/th arc into (side=0) or out of (side=1) /n/. *)
      BEGIN
        <* ASSERT 0 <= i AND i < deg *>
        WITH
          angle = BasicArcDir(i, deg)
        DO
          a.s[1 - side].dir := angle;
        END
      END ProcessArc;

    VAR a: Arc; i: NAT := 0;
    BEGIN
      a := n.s[side].arcs;
      WHILE a # NIL DO
        ProcessArc(a, i);
        INC(i);
        a := a.s[1 - side].next
      END;
    END ProcessNode;
  
  BEGIN
    WITH node = t.node^ DO
      FOR s := 1 TO LAST(node) DO
        WITH n = node[s] DO
          IF n # NIL THEN 
            ProcessNode(n, 0);
            ProcessNode(n, 1)
          END
        END
      END
    END
  END ComputeInitialArcDirs;

PROCEDURE ComputeNodePositionsAndDraw(t: ADDiagram.T; psFile: Wr.T) =

  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    WITH
      node = t.node^,
      gridN = t.gridN
    DO

      PROCEDURE PerformOneIteration(now: NAT): BOOL =
        VAR n: Node;
            range: IntBox;
        BEGIN
          n := ADDiagram.FirstNodeInQueue(t.unplaced);
          IF n = NIL THEN
            n := ADDiagram.FirstNodeInQueue(t.unhappy);
            IF n = NIL THEN
              ReportSuccess();
              RETURN TRUE
            END;
          END;

          IF NOT o.quiet THEN ReportAttemptToPlace(n.label, n.placed) END;
          IF n.placed THEN t.UnplaceNode(n) END;

          range := ComputeCurrentRange(n);
          IF range[0].lo > range[0].hi THEN
            ReportInsufficientHorRange(n.label, range[0]);
            <* ASSERT FALSE *>
          END;
          IF o.verbose THEN ReportCurrentRange(n.label, range) END;

          OptimizeNodePosition(n, range, now := now);
          IF n.gp = Nowhere THEN
            ReportFailureToPlace(n.label);
            <* ASSERT FALSE *>
          ELSE
            t.PlaceNode(n, n.gp, happy := TRUE, now := now)
          END;
          IF NOT o.quiet THEN 
            ReportQueueStatus(QueueLength(t.unplaced), QueueLength(t.unhappy))
          END;
          RETURN FALSE
        END PerformOneIteration;

      PROCEDURE CheckForEnoughSpace() =
        VAR totNodeCells: NAT := 0;  (* Sum of node areas (in cells) *)
        <* FATAL Wr.Failure, Thread.Alerted, BadParameters *>
        BEGIN
          FOR s := 1 TO LAST(node) DO
            IF node[s] # NIL THEN 
              totNodeCells := totNodeCells + ADDiagram.NodeArea(node[s])
            END;
          END;
          IF o.verbose THEN
            Wr.PutText(stderr, "Total node area (cells) = ");
            Wr.PutText(stderr, Fmt.Int(totNodeCells));
            Wr.PutText(stderr, " (average ");
            WITH avg = FLOAT(totNodeCells)/FLOAT(t.nNodes) DO
              Wr.PutText(stderr, Fmt.Real(avg, 2, Fmt.Style.Flo))
            END;
            Wr.PutText(stderr, " per node)\n");
          END;

          WITH
            gridCells = (gridN[0]+1)*(gridN[1]+1)
          DO
            IF totNodeCells > gridCells DIV 3 THEN
              ErrorNotEnoughGridCells(totNodeCells)
            END
          END
        END CheckForEnoughSpace;

      PROCEDURE ComputeCurrentRange(n: Node): IntBox =
        (*
          Computes the min and max coordinates of node /n/,
          taking into account the current positions of its neighbors,
          placed or unplaced.
          *)
        VAR ri: NAT := n.s[0].rad;
            ro: NAT := n.s[1].rad;
            xmin: INT := MAX(ri, n.range[0].lo);            (* Lower bound for /n.gp[0]/ *)
            xmax: INT := MIN(gridN[0] - ro, n.range[0].hi); (* Upper bound for /n.gp[0]/ *)
            
        PROCEDURE ProcessSucc (succ: Node) = 
          BEGIN
            <* ASSERT succ # NIL *>
            IF succ.gp # Nowhere THEN
              WITH
                riSucc = succ.s[0].rad
              DO
                xmax := MIN(xmax, succ.gp[0] - ro - riSucc);
              END
            END
          END ProcessSucc;

        PROCEDURE ProcessPred (pred: Node) =
          BEGIN
            <* ASSERT pred # NIL *>
            IF pred.gp # Nowhere THEN
              WITH
                roPred = pred.s[1].rad
              DO
                xmin := MAX(xmin, pred.gp[0] + ri + roPred)
              END
            END
          END ProcessPred;

        VAR a: Arc;
        BEGIN
          (* Process sucessors: *)
          a := n.s[1].arcs;
          WHILE a # NIL DO
            ProcessSucc(a.s[1].node);
            a := a.s[0].next
          END;

          (* Process predecessors: *)
          a := n.s[0].arcs;
          WHILE a # NIL DO
            ProcessPred(a.s[0].node);
            a := a.s[1].next
          END;

          WITH
            hr = IntRange{xmin, xmax}
          DO
            <* ASSERT hr.lo >= ri *>
            <* ASSERT hr.hi <= gridN[0] - ro *>
            IF n.gp # Nowhere THEN
              <* ASSERT n.gp[0] >= hr.lo AND n.gp[0] <= hr.hi *>
            END;
            RETURN IntBox{hr, n.range[1]}
          END
        END ComputeCurrentRange;

      PROCEDURE CheckPositionConsistency() =
        BEGIN
          FOR s := 1 TO LAST(node) DO
            IF node[s] # NIL THEN 
              WITH range = ComputeCurrentRange(node[s]) DO
                IF o.verbose THEN ReportCurrentRange(node[s].label, range) END;
              END
            END;
          END;
        END CheckPositionConsistency;

      <*UNUSED*>
      PROCEDURE UnplaceNearestNeighbors(n: Node; horRange: IntRange) =
        (*
          Given a node /n/ whose horizontal range /horRange/ is too
          narrow (or empty), attempts to unplace the placed predecessors
          and successors of /n/ that determine that horizontal range, unless
          they are fixed.
          Bombs out if the unplacements done are not enough to widen
          the horizontal range.
          *)

        VAR ri: NAT := n.s[0].rad;
            ro: NAT := n.s[1].rad;
            loPinned: BOOL := (horRange.lo <= MAX(ri, n.range[0].lo));
            hiPinned: BOOL := (horRange.hi >= MIN(gridN[0] - ro, n.range[0].hi));
            
        PROCEDURE ProcessSucc (succ: Node) = 
          BEGIN
            <* ASSERT succ # NIL *>
            IF succ.placed THEN
              <* ASSERT succ.gp[0] # Nowhere[0] *>
              <* ASSERT succ.gp[1] # Nowhere[1] *>
              WITH
                riSucc = succ.s[0].rad
              DO
                IF horRange.hi >= succ.gp[0] - ro - riSucc THEN
                  IF succ.gp[0] >= succ.range[0].hi THEN
                    hiPinned := TRUE
                  ELSE
                    t.UnplaceNode(succ)
                  END
                END
              END
            END
          END ProcessSucc;

        PROCEDURE ProcessPred (pred: Node) =
          BEGIN
            <* ASSERT pred # NIL *>
            IF pred.placed THEN
              <* ASSERT pred.gp[0] # Nowhere[0] *>
              <* ASSERT pred.gp[1] # Nowhere[1] *>
              WITH
                roPred = pred.s[1].rad
              DO
                IF horRange.lo <= pred.gp[0] + ri + roPred THEN
                  IF pred.gp[0] <= pred.range[0].lo THEN
                    loPinned := TRUE
                  ELSE
                    t.UnplaceNode(pred)
                  END
                END
              END
            END
          END ProcessPred;

        VAR a: Arc;
        BEGIN
          (* Process sucessors: *)
          a := n.s[1].arcs;
          WHILE a # NIL DO
            ProcessSucc(a.s[1].node);
            a := a.s[0].next
          END;

          (* Process predecessors: *)
          a := n.s[0].arcs;
          WHILE a # NIL DO
            ProcessPred(a.s[0].node);
            a := a.s[1].next
          END;

          IF hiPinned AND loPinned THEN
            ErrorNotEnoughSpace(n.label)
          END
        END UnplaceNearestNeighbors;

      PROCEDURE OptimizeNodePosition(n: Node; range: IntBox; now: NAT) =
        (*
          Places node "n" at its best position inside "range",
          taking into acount the position of other nodes.

          If node "n" moves, or its arc directions change, its neighbors become 
          unhappy. *)

        VAR
          gpOld, gpIdeal, gpBest: IntPoint;
          cpBest, ppBest: REAL := LAST(REAL);
          minPenalty: REAL := -1.0;

        PROCEDURE TryPosition(gpRaw: IntPoint; kind: TEXT := ""): BOOL =
          (* Assumes arc directions have been recomputed assuming "n" at "n.gp". *)
          BEGIN
            WITH gpNew = ClipIntPoint(gpRaw, range) DO
              IF gpNew # n.gp THEN
                t.MoveNode(n, gpNew);
                RecomputeArcDirections(t, n, saddenNeighbors := FALSE);
              END;
              WITH
                pp = PositionPenalty(t, n, cutoff := ppBest + cpBest),
                cp = CollisionPenalty(t, n, now, cutoff := ppBest + cpBest - pp)
              DO
                IF o.verbose THEN
                  ReportTentativePlacement(n.label, n.gp, cp, pp, kind)
                END;
                IF pp < minPenalty THEN
                  minPenalty := 0.0
                ELSIF minPenalty = -1.0 AND gpNew = gpIdeal THEN
                  minPenalty := pp
                END;
                IF pp + cp < ppBest + cpBest THEN
                  gpBest := n.gp;
                  ppBest := pp;
                  cpBest := cp;
                  RETURN TRUE
                ELSE
                  RETURN FALSE
                END
              END;
            END;
          END TryPosition;
          
        PROCEDURE DecidePosition(gpNew: IntPoint) =
          BEGIN
            t.MoveNode(n, gpNew);
            IF gpNew # gpOld THEN MakeNeighborsUnhappy(n) END;
          END DecidePosition;

        <* FATAL Wr.Failure, Thread.Alerted *>
        BEGIN
          <* ASSERT range[0].lo <= range[0].hi *>
          <* ASSERT range[1].lo <= range[1].hi *>
          
          WITH
            ri = n.s[0].rad,
            ro = n.s[1].rad,
            ht = n.ht
          DO
            gpOld := n.gp;
            IF gpOld # Nowhere THEN
              RecomputeArcDirections(t, n, saddenNeighbors := TRUE);
              IF range[0].lo = range[0].hi AND range[1].lo = range[1].hi THEN
                (* Node can't move, nothing else to do: *)
                DecidePosition(IntPoint{range[0].lo, range[1].lo}); 
                RETURN
              END;
              EVAL TryPosition(gpOld, "current");
            END;
            
            gpIdeal := ComputeIdealPosition(t, n, range);
            IF o.verbose THEN ReportIdealPosition(n.label, gpIdeal) END;
            IF gpIdeal # gpOld THEN
              EVAL TryPosition(gpIdeal, "closest to ideal");
              IF cpBest = 0.0 THEN 
                DecidePosition(gpBest); RETURN
              END;
            END;
            
            CONST
              MaxTrialsNormal = 10;     (* Placements to try, once no collisions. *)
              MaxTrialsDesperate = 50;  (* Placements to try, while there are collisions. *)
            VAR
              maxSize: Counts;   (* Max size of target region where to throw darts *)
              size: Counts;      (* Current size of region where to throw darts *)
              gp: IntPoint;
              trials: NAT := 1;
            BEGIN
              (* maxSize[0] := 1 + ri + ro + 1; *)
              (* maxSize[1] := 1 + 2*ht + 1; *)
              maxSize[0] := range[0].hi - range[0].lo + 1;
              maxSize[1] := range[1].hi - range[1].lo + 1;
              size[0] := 1;
              size[1] := 1;
                            
              WHILE (cpBest > 0.0 AND trials < MaxTrialsDesperate)
              OR (cpBest = 0.0 AND trials < MaxTrialsNormal) DO
                FOR c := 0 TO 1 DO
                  WITH 
                    lo = MAX(gpBest[c] - size[c], range[c].lo),
                    hi = MIN(gpBest[c] + size[c], range[c].hi),
                    t = lo + Random.Subrange(NIL, 0, hi - lo + 1)
                  DO
                    gp[c] := t
                  END;
                END;
                IF TryPosition(gp) THEN
                  size[0] := MAX(1, (2 * size[0] + 2) DIV 3);
                  size[1] := MAX(1, (2 * size[1] + 2) DIV 3);
                ELSE
                  size[0] := MIN(maxSize[0], (3 * size[0] + 1) DIV 2);
                  size[1] := MIN(maxSize[1], (3 * size[1] + 1) DIV 2);
                END;
                INC(trials);
                IF cpBest + ppBest <= minPenalty THEN 
                  DecidePosition(gpBest); RETURN 
                END;
              END;
            END;
            DecidePosition(gpBest)
          END;
        END OptimizeNodePosition;

      PROCEDURE MakeNeighborsUnhappy(n: Node) =
        VAR a: Arc;
        BEGIN
          FOR side := 0 TO 1 DO
            a := n.s[side].arcs;
            WHILE a # NIL DO
              WITH m = a.s[side].node DO
                <* ASSERT m # NIL *>
                IF m.placed AND m.happy THEN t.MakeNodeUnhappy(m) END
              END;
              a := a.s[1-side].next
            END;
          END;
        END MakeNeighborsUnhappy;

      BEGIN
        IF o.verbose THEN PrintGridData(o.gridDims, o.gridStep, gridN) END;
        CheckForEnoughSpace();
        ComputeRanges(t);
        ComputeInitialPositions(t, now := 0);
        CheckPositionConsistency();
        ComputeInitialArcDirs(t);
        
        VAR iterations: NAT := 0;
            success: BOOL := FALSE;
            npages: NAT := 0;
        BEGIN
          WHILE NOT success AND iterations < o.iterations DO
            IF iterations >= o.showStart
            AND npages < o.showCount
            AND (iterations - o.showStart) MOD o.showStep = 0 THEN 
              DrawDiagram(psFile, t, page := npages + 1, iterations := iterations);
              INC(npages);
            END;
            INC(iterations);
            success := PerformOneIteration(now := iterations);
          END;
          DrawDiagram(psFile, t, page := npages + 1, iterations := iterations);
          INC(npages);
          ClosePSFile(psFile, npages := npages);
        END
      END
    END
  END ComputeNodePositionsAndDraw;
  
PROCEDURE ComputeRanges(t: ADDiagram.T) =
  (*
    Computes the maximum range "n.range" for the position of every node "n",
    allowing for any rearrangement of all nodes.
    Takes into account only the grid limits and node sizes. *)
  
  PROCEDURE ComputeMaxHorPosition(n: Node) =
    VAR ro: NAT := n.s[1].rad;
        maxH: INT := t.gridN[0] - ro;
        a: Arc;
    BEGIN
      a := n.s[1].arcs;
      WHILE a # NIL DO
        WITH succ = a.s[1].node DO
          <* ASSERT succ # NIL *>
          WITH
            riSucc = succ.s[0].rad,
            maxHSucc = succ.range[0].hi
          DO
            <* ASSERT maxHSucc <= t.gridN[0] *>
            maxH := MIN(maxH, maxHSucc - riSucc - ro);
          END
        END;
        a := a.s[0].next
      END;
      n.range[0].hi := maxH
    END ComputeMaxHorPosition;

  PROCEDURE ComputeMinHorPosition(n: Node) =
    VAR ri: NAT := n.s[0].rad;
        minH: INT := ri;
        a: Arc;
    BEGIN
      a := n.s[0].arcs;
      WHILE a # NIL DO
        WITH pred = a.s[0].node DO
          <* ASSERT pred # NIL *>
          WITH
            roPred = pred.s[1].rad,
            minHPred = pred.range[0].lo
          DO
            <* ASSERT minHPred >= 0 *>
            minH := MAX(minH, minHPred + roPred + ri);
          END
        END;
        a := a.s[1].next
      END;
      n.range[0].lo := minH
    END ComputeMinHorPosition;

  <* FATAL Wr.Failure, Thread.Alerted, BadParameters *>
  BEGIN
    WITH
      node = t.node^
    DO
      (* Computes maximum x position for every node: *)
      FOR s := 1 TO LAST(node) DO 
        WITH n = node[s] DO 
          IF n # NIL THEN ComputeMaxHorPosition(n) END
        END
      END;

      (* Computes minimum x position for every node: *)
      FOR s := LAST(node) TO 1 BY -1 DO 
        WITH n = node[s] DO 
          IF n # NIL THEN ComputeMinHorPosition(n) END
        END
      END;
      
      (* Sets vertical range for every node: *)
      FOR s := 1 TO LAST(node) DO 
        WITH n = node[s] DO 
          IF n # NIL THEN
            n.range[1] := IntRange{n.ht, t.gridN[1] - n.ht};
          END;
        END
      END;
      
      (* Checks if constraints are satisfiable: *)
      FOR s := 1 TO LAST(node) DO 
        WITH n = node[s] DO 
          IF n # NIL THEN
            IF n.range[0].lo > n.range[0].hi 
            OR n.range[1].lo > n.range[1].hi 
            THEN 
              ErrorNotEnoughSpace(n.label)
            END;
          END;
        END
      END;
    END
  END ComputeRanges;

PROCEDURE XYFromGP(i: INT): REAL =
  BEGIN
    RETURN FLOAT(i) * o.gridStep
  END XYFromGP;

<*UNUSED*>
PROCEDURE GPFromXY(READONLY x: REAL): INT =
  BEGIN
    RETURN ROUND(x / o.gridStep)
  END GPFromXY;

PROCEDURE NaturalArcDirection(a: Arc; side: [0..1]): REAL =
  (* 
    Computes the direction of arc "a" at its "side" endpoint, 
    as if it were the only arc incident to the node at that end *)
  BEGIN
    WITH 
      org = a.s[0].node,
      dst = a.s[1].node 
    DO
      <* ASSERT org # NIL *>
      <* ASSERT dst # NIL *>
      IF org.gp[0] = Nowhere[0] OR org.gp[1] = Nowhere[1]
      OR dst.gp[0] = Nowhere[0] OR dst.gp[1] = Nowhere[1]
      THEN
        RETURN 0.0
      ELSE
        WITH
          dx = FLOAT(dst.gp[0] - org.gp[0]),
          dy = FLOAT(dst.gp[1] - org.gp[1]),
          dirStraight = Dir(dx, dy),
          dirBent = dirStraight - 0.5 * (a.s[1-side].dir - dirStraight)
        DO
          IF ADDiagram.Degree(a.s[1-side].node, side) = 1 THEN
            RETURN dirStraight
          ELSE
            RETURN dirBent
          END
        END
      END;
    END
  END NaturalArcDirection;

PROCEDURE RecomputeArcDirections(t: ADDiagram.T; n: Node; saddenNeighbors: BOOL) =
  (*
    Recomputes the nominal arriving and departing directions (".dir" fields)
    of all arcs incident to node "n", and also "n.tilt".
    
    For every arc "a", if "saddenNeighbors" is true and "n.tilt" 
    or "a.dir" have changed, makes the other endpoint of "a" unhappy. *)

  VAR degs: ARRAY [0..1] OF NAT;
      arcs: ARRAY [0..1] OF REF ARRAY OF Arc;
      dirs: ARRAY [0..1] OF REF ARRAY OF REAL;

  PROCEDURE AllocVariables(side: [0..1]) =
    BEGIN
      degs[side] := ADDiagram.Degree(n, side);
      arcs[side] := NEW(REF ARRAY OF Arc, degs[side]);
      dirs[side] := NEW(REF ARRAY OF REAL, degs[side]);
    END AllocVariables;

  PROCEDURE ComputeRawSideDirs(side: [0..1]) =
    BEGIN
      WITH
        deg = degs[side],
        arc = arcs[side]^,
        dir = dirs[side]^
      DO
        (* Collect arcs: *)
        VAR a: Arc := n.s[side].arcs; i: NAT := 0;
        BEGIN
          WHILE a # NIL DO 
            arc[i] := a; 
            dir[i] := NaturalArcDirection(a, 1-side);
            INC(i); 
            a := a.s[1-side].next 
          END;
          <* ASSERT i = deg *>
        END;
      END
    END ComputeRawSideDirs;
    
  PROCEDURE SortArcsByDir(side: [0..1]) =
    BEGIN
      WITH
        deg = degs[side],
        arc = arcs[side]^,
        dir = dirs[side]^
      DO
        IF deg > 1 THEN
          (* Sort all arcs by increasing direction: *)
          FOR i := 1 TO deg-1 DO
            VAR j := i;
            BEGIN
              WHILE j > 0 AND dir[j-1] > dir[j] DO
                VAR ta: Arc; td: REAL; 
                BEGIN 
                  ta := arc[j-1]; arc[j-1] := arc[j]; arc[j] := ta; 
                  td := dir[j-1]; dir[j-1] := dir[j]; dir[j] := td; 
                END;
                DEC(j);
              END
            END
          END;
        END;
      END
    END SortArcsByDir;
    
  PROCEDURE ComputeMeanDir(): REAL =
    VAR sum: REAL := 0.0;
    BEGIN
      IF degs[0] = 0 AND degs[1] = 0 THEN
        RETURN 0.0
      ELSE
        FOR side := 0 TO 1 DO
          FOR i := 0 TO degs[side]-1 DO 
            sum := sum + dirs[side][i]
          END;
        END;
        RETURN sum / FLOAT(degs[0] + degs[1])
      END
    END ComputeMeanDir;

  PROCEDURE AdjustSideDirs(side: [0..1]) =
    (* Recomputes "dir" table based on rank: *)
    BEGIN
      IF degs[side] # 0 THEN
        WITH
          deg = degs[side],
          arc = arcs[side]^,
          dir = dirs[side]^
        DO
          FOR i := 0 TO deg-1 DO 
            WITH
              frac = FLOAT(i+1)/FLOAT(deg+1),
              d = (frac - 0.5) * MaxFanWidth
            DO
              dir[i] := d
            END
          END;
        END
      END
    END AdjustSideDirs;
    
  PROCEDURE ComputeTilt(meanDir: REAL): REAL =
    CONST TiltGrain = 0.5*Degree;
    BEGIN
      IF degs[0] = 0 OR degs[1] = 0 THEN
        RETURN 0.0
      ELSE
        WITH
          minDir = MIN(dirs[0][0], dirs[1][0]),
          maxDir = MAX(dirs[0][degs[0]-1], dirs[1][degs[1]-1]),
          minTilt = - 0.5*MaxFanWidth - minDir,
          maxTilt = + 0.5*MaxFanWidth - maxDir,
          rawTilt = MAX(minTilt, MIN(maxTilt, meanDir)),
          tilt = TiltGrain * FLOAT(TRUNC(rawTilt/TiltGrain))
        DO
          <* ASSERT minTilt <= 0.0 *>
          <* ASSERT maxTilt >= 0.0 *>
          RETURN tilt
        END
      END
    END ComputeTilt;
    
  PROCEDURE SetSideDirs(side: [0..1]) =
    BEGIN
      WITH
        deg = degs[side],
        arc = arcs[side]^,
        dir = dirs[side]^
      DO
        (* Set ".dir" fields: *)
        FOR i := 0 TO deg-1 DO
          WITH a = arc[i], d = dir[i] DO 
            IF a.s[1-side].dir # d THEN
              a.s[1-side].dir := d;
              WITH m = a.s[side].node DO
                IF m.placed AND saddenNeighbors THEN t.MakeNodeUnhappy(m) END
              END
            END;
          END;
        END
      END;
    END SetSideDirs;
    
  PROCEDURE SetTilt(newTilt: REAL) =
    (* Set "n.tilt": *)
    BEGIN
      IF n.tilt # newTilt THEN
        IF saddenNeighbors THEN 
          FOR side := 0 TO 1 DO
            WITH
              deg = degs[side],
              arc = arcs[side]^
            DO
              FOR i := 0 TO deg-1 DO
                WITH a = arc[i], m = a.s[side].node DO
                  IF m.placed AND m.happy THEN t.MakeNodeUnhappy(m) END
                END
              END
            END;
          END;
        END;
        n.tilt := newTilt;
      END;
    END SetTilt;
  
  VAR meanDir, newTilt: REAL;
  BEGIN
    FOR side := 0 TO 1 DO 
      AllocVariables(side);
      ComputeRawSideDirs(side);
    END;
    meanDir := ComputeMeanDir();
    FOR side := 0 TO 1 DO 
      SortArcsByDir(side);
      AdjustSideDirs(side)
    END;
    newTilt := ComputeTilt(meanDir);
    FOR side := 0 TO 1 DO 
      SetSideDirs(side)
    END;
    SetTilt(newTilt);
  END RecomputeArcDirections;

PROCEDURE ComputeArcLabelCoords(a: Arc; READONLY cp: ArcControlPoints): RealPoint =
  BEGIN
    WITH
      org = a.s[0].node,
      
      dir0 = a.s[0].dir,
      
      xa = cp[1][0] * o.gridStep,
      ya = cp[1][1] * o.gridStep,

      xlab = xa - 0.5 * o.arcSpacing * Sin(dir0),
      ylab = ya + 0.5 * o.arcSpacing * Cos(dir0)
    DO
      RETURN RealPoint{xlab, ylab}
    END
  END ComputeArcLabelCoords;

PROCEDURE OpenPSFile (): Wr.T =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    WITH
      wr = FileStream.OpenWrite(o.autName & ".ps")
    DO
      ADPlot.BeginFile(wr);
      RETURN wr
    END
  END OpenPSFile;

PROCEDURE ClosePSFile (wr: Wr.T; npages: NAT) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    ADPlot.EndFile(wr, npages);
    Wr.Close(wr);
  END ClosePSFile;

PROCEDURE DrawDiagram(wr: Wr.T; t: ADDiagram.T; page: NAT := 99; iterations: NAT := LAST(NAT)) =

  PROCEDURE DrawNode(n: Node) =
    <* FATAL Wr.Failure, Thread.Alerted *>
    VAR gray: REAL := 1.0;
    BEGIN
      IF n.gp[0] # Nowhere[0] AND n.gp[1] # Nowhere[1] THEN
        WITH 
          x = L(XYFromGP(n.gp[0])),
          y = L(XYFromGP(n.gp[1])),
          r = L(o.nodeRadius)
        DO
          IF NOT n.placed THEN gray := 0.75 END;
          ADPlot.FillAndDrawCircle(wr, x, y, r, gray);
          IF n.marked THEN
            ADPlot.FillAndDrawCircle(wr, x, y, 0.8D0*r, gray)
          END;
          ADPlot.PutLabel(wr, Fmt.Int(n.label), x, y, 0.5, 0.5);
        END;
      END;
    END DrawNode;

  PROCEDURE DrawArc (a: Arc) =
    <* FATAL Wr.Failure, Thread.Alerted *>
    BEGIN
      WITH
        org = a.s[0].node,
        dst = a.s[1].node
      DO      
        IF org.gp[0] # Nowhere[0] OR org.gp[1] # Nowhere[1]
        OR dst.gp[0] # Nowhere[0] OR dst.gp[1] # Nowhere[1] THEN
          WITH
            cp = ADDiagram.ComputeArcControlPoints(a)
          DO
            (* Draw curve: *)
            WITH
              x0 = L(cp[0][0] * o.gridStep),
              y0 = L(cp[0][1] * o.gridStep),

              xa = L(cp[1][0] * o.gridStep),
              ya = L(cp[1][1] * o.gridStep),

              xb = L(cp[2][0] * o.gridStep),
              yb = L(cp[2][1] * o.gridStep),

              xc = L(cp[3][0] * o.gridStep),
              yc = L(cp[3][1] * o.gridStep),

              xd = L(cp[4][0] * o.gridStep),
              yd = L(cp[4][1] * o.gridStep)
            DO
              ADPlot.DrawSegment(wr, x0, y0, xa, ya);
              ADPlot.DrawCurve(wr, xa, ya, xb, yb, xc, yc, xd, yd);
            END;
            (* Print label: *)
            WITH
              labp = ComputeArcLabelCoords(a, cp),
              labt = Text.FromChar(a.label)
            DO
              ADPlot.PutLabel(wr, labt, L(labp[0]), L(labp[1]), 0.5, 0.5)
            END;
          END                
        END
      END;
    END DrawArc;

  <* FATAL ADDiagram.Abort *>
  BEGIN

    ADPlot.BeginPage(wr, page, 1, 1);
    ADPlot.AddCaption(wr, "automaton " & o.autName);
    ADPlot.AddCaption(wr, "pass " & Fmt.Int(iterations));

    (* Draw arcs: *)
    ADPlot.BeginSection(wr, "Arcs");
    ADPlot.SetLabelFontSize(wr, 5.0);
    t.EnumArcs(DrawArc);
    ADPlot.EndSection(wr);

    (* Draw nodes: *)
    ADPlot.BeginSection(wr, "Nodes");
    ADPlot.SetLabelFontSize(wr, 5.0);
    t.EnumNodes(DrawNode);
    ADPlot.EndSection(wr);

    ADPlot.DrawFrame(wr);
    ADPlot.EndPage(wr);

  END DrawDiagram;

PROCEDURE FT8(x: REAL): TEXT =
  BEGIN
    RETURN Fmt.Pad(Fmt.Real(x, 2, Fmt.Style.Flo), 8)
  END FT8;

PROCEDURE WriteCoordFiles(t: ADDiagram.T) =
  <* FATAL Wr.Failure, Rd.Failure, Thread.Alerted *>
  BEGIN
    WITH
      wr = FileStream.OpenWrite(o.autName & ".nxy")
    DO

      PROCEDURE PrintNode(n: Node) =
        <* FATAL Wr.Failure, Thread.Alerted *>
        BEGIN
          WITH
            x = XYFromGP(n.gp[0]),
            y = XYFromGP(n.gp[1])
          DO
            Wr.PutText(wr, Fmt.Pad(Fmt.Int(n.label), 14));
            Wr.PutText(wr, "  ");
            Wr.PutText(wr, FT8(x));
            Wr.PutText(wr, " ");
            Wr.PutText(wr, FT8(y));
            Wr.PutText(wr, " \n");
          END
        END PrintNode;
        
      <* FATAL ADDiagram.Abort *>
      BEGIN
        t.EnumNodes(PrintNode)
      END;
      Wr.Close(wr)
    END;

    WITH
      wr = FileStream.OpenWrite(o.autName & ".axy")
    DO
    
      PROCEDURE PrintArc(a: Arc) =
        <* FATAL Wr.Failure, Thread.Alerted *>
        BEGIN
          WITH
            cp = ADDiagram.ComputeArcControlPoints(a),
            lp = ComputeArcLabelCoords(a, cp)
          DO
            Wr.PutText(wr, "(");
            Wr.PutChar(wr, a.label);
            Wr.PutText(wr, ")");
            Wr.PutText(wr, " ");
            Wr.PutText(wr, FT8(lp[0]));
            Wr.PutText(wr, " ");
            Wr.PutText(wr, FT8(lp[1]));
            Wr.PutText(wr, " ");
            FOR i := 0 TO LAST(cp) DO
              Wr.PutText(wr, "  ");
              Wr.PutText(wr, FT8(cp[i][0] * o.gridStep));
              Wr.PutText(wr, " ");
              Wr.PutText(wr, FT8(cp[i][1] * o.gridStep));
            END;
            Wr.PutText(wr, "\n");
          END
        END PrintArc;

      <* FATAL ADDiagram.Abort *>
      BEGIN
        t.EnumArcs(PrintArc);
      END;
      Wr.Close(wr)
    END;
  END WriteCoordFiles;
  
EXCEPTION
  BadParameters;    (* Bad input parameters *)

PROCEDURE GetOptions(): Options =
  VAR o: Options;
  <* FATAL Scan.BadFormat, Wr.Failure, Thread.Alerted *>
  BEGIN
    ParseParams.BeginParsing(stderr);

      ParseParams.GetKeyword("-autName");
      o.autName := ParseParams.GetNext();

      o.writeCoords := ParseParams.KeywordPresent("-writeCoords");

      o.verbose := ParseParams.KeywordPresent("-verbose");
      o.quiet :=  ParseParams.KeywordPresent("-quiet");

      IF ParseParams.KeywordPresent("-gridDims") THEN
        o.gridDims[0] := ParseParams.GetNextReal(0.01, 100.0);
        o.gridDims[1] := ParseParams.GetNextReal(0.01, 100.0);
      ELSE
        o.gridDims[0] := FLOAT(ADPlot.XMax, REAL);
        o.gridDims[1] := FLOAT(ADPlot.YMax, REAL);
      END;

      IF ParseParams.KeywordPresent("-iterations") THEN
        o.iterations := ParseParams.GetNextInt(1, LAST(NAT));
      ELSE
        o.iterations := 1000;
      END;

      IF ParseParams.KeywordPresent("-showStart") THEN
        o.showStart := ParseParams.GetNextInt(0, LAST(NAT));
      ELSE
        o.showStart := 0;
      END;

      IF ParseParams.KeywordPresent("-showCount") THEN
        o.showCount := ParseParams.GetNextInt(1, 1000);
      ELSE
        o.showCount := 50;
      END;

      IF ParseParams.KeywordPresent("-showStep") THEN
        o.showStep := ParseParams.GetNextInt(1, LAST(NAT));
      ELSE
        o.showStep := ((o.iterations - o.showStart) + o.showCount - 1) DIV o.showCount;
      END;

      IF ParseParams.KeywordPresent("-arcSpacing") THEN
        o.arcSpacing := ParseParams.GetNextReal(0.001, 0.010);
      ELSE
        o.arcSpacing := 0.003;
      END;

      IF ParseParams.KeywordPresent("-nodeRadius") THEN
        o.nodeRadius := ParseParams.GetNextReal(0.001, 0.010);
      ELSE
        o.nodeRadius := 0.003;
      END;

      WITH md = MAX(o.gridDims[0], o.gridDims[1]) DO
        IF ParseParams.KeywordPresent("-gridStep") THEN
          o.gridStep := ParseParams.GetNextReal(
            md/FLOAT(ADDiagram.MaxGridN), 
            md/FLOAT(ADDiagram.MinGridN)
          );
        ELSE
          o.gridStep := md/FLOAT(ADDiagram.MaxGridN)
        END;
      END;

    ParseParams.EndParsing();
    RETURN o
  END GetOptions;

BEGIN
  Main()
END AutoDrawMain.

