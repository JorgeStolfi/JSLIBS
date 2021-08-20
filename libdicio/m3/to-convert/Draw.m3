(****************************************************************************)
(* (C) Copyright 1992 Universidade Estadual de Campinas (UNICAMP)           *)
(*                    Campinas, SP, Brazil                                  *)
(*                                                                          *)
(* Authors:                                                                 *)
(*                                                                          *)
(*   Tomasz Kowaltowski  - CS Dept, UNICAMP <tomasz@dcc.unicamp.ansp.br>    *)
(*   Claudio L. Lucchesi - CS Dept, UNICAMP <lucchesi@dcc.unicamp.ansp.br>  *)
(*   Jorge Stolfi        - DEC Systems Research Center <stolfi@src.dec.com> *)
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

MODULE Draw EXPORTS Main;

IMPORT 
  Rd, Wr, Char, Text, Fmt, FileStream, ParseParams, 
  Math, Random,
  Basics, StringPrinter, Reduced, ReducedPair;
FROM Basics IMPORT INT, NAT, POS, BOOL;
FROM Basics IMPORT Done, Skip, Abort, WrNat;
FROM Reduced IMPORT Letter, String, State, Arc, NullState;
FROM ReducedPair IMPORT BoolOp, Which;
FROM Stdio IMPORT stdin, stdout, stderr;

CONST
  MaxGridN = 2000;       (* Maximum number of positions in grid, each axis *)
  MinGridN = 20;         (* Minimum number of positions in grid, each axis *)
  MaxNodes = 256*256-1;  (* Maximum number of nodes that can be placed *)
  
TYPE
  Node = [0..MaxNodes];  (* Number of a proper node, or NullNode (0) *)

CONST
  NullNode = 0;          (* A dummy node number *)
  
CONST
  Nowhere = IntPoint{LAST(NAT), LAST(NAT)};

TYPE
  Drawing = REF RECORD
      nodeData: REF ARRAY OF NodeData; (* Maps Node numbers to node data *)
      node: REF ARRAY OF Node;         (* Maps State numbers to Node numbers *)
      gDims: Dims;                     (* Grid dimensions (meters) *)
      gStep: REAL;                     (* Grid mesh size (meters) *)
    END;

  NodeData = RECORD
      state: State;    (* State number in automaton *)
      gp: IntPoint;    (* Position in grid *)
      ni: NAT;         (* In-degree or state *)
      ri: REAL;        (* Radius of input territory *)
      no: NAT;         (* Out-degree or state *)
      ro: REAL;        (* Radius of output territory *)
    END;

TYPE

  IntPoint = ARRAY [0..1] OF INT;     (* Indices of grid point *)

  RealPoint = ARRAY [0..1] OF REAL;   (* Cartesian coordinates *)
  
  RealRange = RECORD lo, hi: REAL END;  (* Set OF REALs between /lo/ and /hi/, incl. *)
  
  IntRange = RECORD lo, hi: INT END;    (* Set of INTs between /lo/ and /hi/, incl. *)

  Dims = ARRAY [0..1] OF REAL;       (* Pair of dimensions *)
  
  Counts = ARRAY [0..1] OF NAT;      (* Pair of counts *)

TYPE
  Options = RECORD
      loadFileName: TEXT;   (* Name of automaton dump-file to read *)
      printFileName: TEXT;  (* Name of file for node positions, or "" if none. *)
      psFileName: TEXT;     (* Name of output Postscript file, or "" if none. *)
      verbose: BOOL;        (* TRUE to make noises while thinking. *)
      gridDims: Dims;       (* Size of drawing area *)
      gridStep: REAL;       (* Mesh of grid of allowed state positions *)
      arcSpacing: REAL;     (* Spacing between arcs at the output radius *)
    END;

TYPE 
  PackedNode = BITS 16 FOR Node;
  PackedBool = BITS 8 FOR BOOL;

VAR o: Options;

PROCEDURE Main() =

  VAR aut: Reduced.T;
      
  BEGIN
    o := GetOptions();
    
    (* Load the automaton: *)
    <* ASSERT NOT Text.IsEmpty(o.loadFileName) *>
    Wr.PutText(stderr, "\n=== loading automaton from " & o.loadFileName & " ===\n\n");
    WITH rd = FileStream.OpenRead(o.loadFileName) DO
      aut := Reduced.Load(rd)
    END;

    (* Print report *)
    WITH ct = aut.Count(ARRAY OF State{aut.Root()}) DO
      Wr.PutText(stderr, "root state = " & Fmt.Pad(Fmt.Int(aut.Root()), 8) & "\n");
      Wr.PutText(stderr, "states =     " & Fmt.Pad(Fmt.Int(ct.states), 8) & "\n");
      Wr.PutText(stderr, "finals =     " & Fmt.Pad(Fmt.Int(ct.finals), 8) & "\n");
      Wr.PutText(stderr, "arcs =       " & Fmt.Pad(Fmt.Int(ct.arcs), 8) & "\n");
    END;
    
    Wr.PutText(stderr, "o.gridStep = " & Fmt.Real(o.gridStep) & "\n");

    WITH
      d = DrawAutomaton(
        aut := aut, 
        state := CollectReachableStates(aut)^,
        gDims := o.gridDims, 
        gStep := o.gridStep
      )
    DO
      IF NOT Text.IsEmpty(o.printFileName) THEN
        WITH
          wr = FileStream.OpenWrite(o.printFileName)
        DO
          PrintDrawing(wr, d);
          Wr.Flush(wr);
        END;
      END;
    END;
    
    Wr.Flush(stderr);
    Wr.Flush(stdout);
  END Main;
  
PROCEDURE PrintDrawing(wr: Wr.T; d: Drawing) =
  BEGIN
    WITH 
      nodeData = d.nodeData^,
      node = d.node^,
      gDims = d.gDims,
      gStep = d.gStep
    DO
      PROCEDURE PrintAllNodes() =
        BEGIN
          FOR n := 1 TO LAST(nodeData) DO
            PrintNode(n);
          END;
        END PrintAllNodes;
        
      PROCEDURE PrintNode(n: Node) =
        BEGIN
          Wr.PutText(wr, "node ");
          Wr.PutText(wr, Fmt.Pad(Fmt.Int(n), 6));
          Wr.PutText(wr, " = state ");
          Wr.PutText(wr, Fmt.Pad(Fmt.Int(nodeData[n].state), 6));
          Wr.PutText(wr, " at ");
          PrintIntPoint(wr, nodeData[n].gp);
          Wr.PutText(wr, "\n");
        END PrintNode;
        
      PROCEDURE PrintAllArcs() =
        BEGIN
          FOR n := 1 TO LAST(nodeData) DO
            PrintNodeArcs(n);
          END;
        END PrintAllArcs;
        
      PROCEDURE PrintNodeArcs(<*UNUSED*> n: Node) =
        BEGIN
        END PrintNodeArcs;
        
      BEGIN
        PrintAllNodes();
        PrintAllArcs();
      END
      
    END
  END PrintDrawing;

PROCEDURE CollectReachableStates(aut: Reduced.T): REF ARRAY OF State =
(*
  Returns a vector with all true states of the automaton /aut/
  reachable from /aut.Root()/, in increasing numerical order. *)
  VAR n: NAT := 0;
  BEGIN
    WITH
      nStates = aut.Count(ARRAY OF State{aut.Root()}).states,
      rres = NEW(REF ARRAY OF State, nStates),
      res = rres^
    DO
      FOR s := 1 TO aut.Root() DO
        IF aut.NPrefs(s) > 0 THEN
          res[n] := s;
          INC(n)
        END
      END;
      <* ASSERT n = nStates *>
      RETURN rres
    END
  END CollectReachableStates;

PROCEDURE DrawAutomaton(
    aut: Reduced.T;                    (* The automaton to draw *)
    READONLY state: ARRAY OF State;    (* States to draw, in increasing order *)
    gDims: Dims;                       (* Grid dimensions *)
    READONLY gStep: REAL;              (* Grid step *)
  ): Drawing = 

  BEGIN
    <* ASSERT NUMBER(state) > 0 *>

    Wr.PutText(stderr, "gStep = " & Fmt.Real(gStep) & "\n");


    WITH
      gridN = Counts{
        ROUND(gDims[0]/gStep), 
        ROUND(gDims[1]/gStep)
      },
      maxNode = NUMBER(state),
      minState = state[0],
      maxState = state[LAST(state)],
      drawing = NEW(Drawing,
        nodeData := NEW(REF ARRAY OF NodeData, maxNode + 1),
        node := NEW(REF ARRAY OF Node, maxState + 1),
        gDims := gDims,
        gStep := gStep
      ),
      node = drawing.node^,
      nodeData = drawing.nodeData^,
      placed = NEW(REF ARRAY OF PackedBool, maxNode + 1)^,
      fixed = NEW(REF ARRAY OF PackedBool, maxNode + 1)^,
      when = NEW(REF ARRAY OF NAT, maxNode + 1)^,
      who = NEW(REF ARRAY OF ARRAY OF PackedNode, gridN[0] + 1, gridN[1] + 1)^,
      nodeQueue = NEW(REF ARRAY OF Node, maxNode + 2)^
    DO
      
      PROCEDURE PlaceNodes(effort: NAT) =
        VAR n: Node;
            gpIdeal, gpBest: IntPoint;
            horRange: RealRange;
        BEGIN
          GatherNodeData();
          CheckForEnoughSpace();
          ClearWhoMap();
          MakeInitialPlacements(now := 0);
          FOR now := 1 TO effort DO
            n := TakeUnplacedNodeFromQueue();
            IF n = NullNode THEN
              ReportSuccess();
              RETURN
            END;
            ReportAttemptToPlace(n);
            
            horRange := ComputeHorRange(n);
            WHILE CEILING(horRange.lo/gStep) > FLOOR(horRange.hi/gStep) DO
              ReportInsufficientHorRange(n, horRange);
              UnplaceNearestNeighbors(n, horRange);
              horRange := ComputeHorRange(n)
            END;

            gpIdeal := IdealPosition(n);
            ReportIdealPosition(n, gpIdeal);
            
            gpBest := ComputeBestPlacement(n, gpIdeal, horRange);
            IF gpBest = Nowhere THEN
              ReportFailureToPlace(n);
              PutUnplacedNodeInQueue(n)
            ELSE
              UnplaceConflictingNodes(n, gpBest);
              ReportPlacement(n, gpBest);
              PlaceNode(n, gpBest, now := now);
            END;
            ReportStatus();
          END
        END PlaceNodes;
        
      PROCEDURE GatherNodeData() =
      (* 
        Gathers/computes the parameters of the nodes to be drawn:
        *)
        BEGIN
          node[NullState] := NullNode;
          nodeData[NullNode] := NodeData{
            state := NullState, 
            gp := Nowhere, 
            ni := 0, ri := 0.0,
            no := 0, ro := 0.0
          };

          IF o.verbose THEN Wr.PutText(stderr, "Collecting node data...\n") END;

          FOR s := 0 TO LAST(node) DO node[s] := NullNode END;
          
          FOR n := 1 TO maxNode DO
            WITH 
              s = state[n-1],
              ni = aut.InDeg(s),
              no = aut.OutDeg(s)
            DO
              node[s] := n;
              nodeData[n] := NodeData{
                state := s,
                gp := Nowhere,
                no := no,
                ro := OutRadius(no),
                ni := ni,
                ri := InRadius(ni)
              };
            END;
          END;
        END GatherNodeData;
        
      PROCEDURE CheckForEnoughSpace() =
        VAR totNodeCells: NAT := 0;  (* Sum of node areas (in cells) *)
        BEGIN
          FOR n := 1 TO maxNode DO 
            totNodeCells := totNodeCells + NodeCells(n)
          END;
          IF o.verbose THEN
            Wr.PutText(stderr, "Total node area (cells) = ");
            Wr.PutText(stderr, Fmt.Int(totNodeCells));
            Wr.PutText(stderr, " (average ");
            WITH avg = FLOAT(totNodeCells)/FLOAT(maxNode) DO
              Wr.PutText(stderr, Fmt.Real(avg, 2, Fmt.Style.Flo))
            END;
            Wr.PutText(stderr, " per node)\n");
          END;
            
          WITH
            gridCells = (gridN[0]+1)*(gridN[1]+1)
          DO
            IF totNodeCells > gridCells DIV 3 THEN
              Wr.PutText(stderr, "** Total node area is too big for given grid **\n");
              RAISE BadParameters
            END
          END
        END CheckForEnoughSpace;
        
      PROCEDURE NodeCells(n: Node): NAT =
      (*
        Area of node /n/, in grid cells (i.e. how many entries of /who/
        it woudl cover if placed).
        *)
        BEGIN
          WITH 
            iri = CEILING(nodeData[n].ri),
            iro = CEILING(nodeData[n].ro),
            iht = MAX(iri, iro)
          DO
            RETURN (iri + 1 + iro) * (iht + 1 + iht)
          END
        END NodeCells;
      
      PROCEDURE MakeInitialPlacements(now: NAT) =
      (*
        Places one "seed" node (possibly more) and puts rest
        in the unplaced node queue. *)
        VAR gp: IntPoint;
            fix: BOOL;
        BEGIN
          FOR n := 1 TO maxNode DO 
            
            (* Decide position /gp/ and whether it is fixed or mobile: *)
            IF nodeData[n].state = aut.Root() THEN
              (* Root node, place at middle of left edge: *)
              gp := IntPoint{0, gridN[1] DIV 2};
              fix := TRUE;
            ELSIF nodeData[n].state = 1 THEN
              (* Unit node, place at middle of right edge: *)
              gp := IntPoint{gridN[0], gridN[1] DIV 2};
              fix := TRUE;
            ELSIF n = maxNode DIV 2 THEN
              (* Middling node, place initially at center of page: *)
              gp := IntPoint{gridN[0] DIV 2, gridN[1] DIV 2};
              fix := FALSE;
            ELSE
              gp := Nowhere;
              fix := FALSE
            END;
            
            (* Now place OR put in unplaced queue, as appropriate: *)
            IF gp # Nowhere THEN
              UnplaceConflictingNodes(n, gp);
              PlaceNode(n, gp, now := now);
              placed[n] := TRUE;
              fixed[n] := fix;
            ELSE
              placed[n] := FALSE;
              fixed[n] := FALSE;
              PutUnplacedNodeInQueue(n)
            END;
          END;
          
        END MakeInitialPlacements;

      CONST MinRadius = 0.003;
      
      PROCEDURE InRadius(<*UNUSED*> indeg: NAT): REAL =
        BEGIN
          RETURN MinRadius;
        END InRadius;
        
      PROCEDURE OutRadius(outdeg: NAT): REAL =
        BEGIN
          RETURN (FLOAT(outdeg + 1) * o.arcSpacing) / 3.1415926
        END OutRadius;

      PROCEDURE IdealPositionFromStateNumber(s: State): IntPoint =
      (*
        Computes an initial grid position for a
        node corresponding to state /s/, based only on the state's number. 
        *)
        VAR x: REAL;
            y: REAL;
        BEGIN
          IF maxState = minState THEN
            x := 0.5 * gDims[0]
          ELSE
            x := gDims[0] * (1.0 - FLOAT(s-minState)/FLOAT(maxState-minState))
          END;
          y := 0.5 * gDims[1];
          RETURN IntPoint{GPFromXY(x), GPFromXY(y)}
        END IdealPositionFromStateNumber;

      PROCEDURE ClearWhoMap() =
      (*
        Clears the /who/ map:
        *)
        BEGIN
          IF o.verbose THEN Wr.PutText(stderr, "Clearing /who/ map...") END;
          FOR x := 0 TO gridN[0] DO
            FOR y := 0 TO gridN[1] DO
              who[x,y] := NullNode
            END
          END;
          IF o.verbose THEN Wr.PutText(stderr, "\n") END;
        END ClearWhoMap;

      VAR frontQ: NAT := 0;
          rearQ: NAT := 0;
      
      PROCEDURE TakeUnplacedNodeFromQueue(): Node =
        BEGIN
          IF frontQ = rearQ THEN 
            RETURN NullNode
          ELSE
            WITH n = nodeQueue[frontQ] DO
              <* ASSERT n # NullNode *>
              INC(frontQ);
              IF frontQ > LAST(nodeQueue) THEN frontQ := 0 END;
              RETURN n
            END
          END
        END TakeUnplacedNodeFromQueue;
      
      PROCEDURE PutUnplacedNodeInQueue(n: Node) =
        BEGIN
          <* ASSERT n # NullNode *>
          nodeQueue[rearQ] := n;
          INC(rearQ);
          IF rearQ > LAST(nodeQueue) THEN rearQ := 0 END;
          <* ASSERT rearQ # frontQ *>
        END PutUnplacedNodeInQueue;
      
      PROCEDURE ReportStatus() =
        BEGIN
          WITH
            sizeQ = NUMBER(nodeQueue),
            nQ = (rearQ + sizeQ - frontQ) MOD sizeQ
          DO
            Wr.PutText(stderr, "  unplaced nodes = " & Fmt.Int(nQ) & "\n")
          END
        END ReportStatus;

      PROCEDURE ReportSuccess() =
        BEGIN
          Wr.PutText(stderr, "Success!\n")
        END ReportSuccess;
        
      PROCEDURE ReportAttemptToPlace(n: Node) =
        BEGIN
          Wr.PutChar(stderr, '\n');
          Wr.PutText(stderr, "attempting to place node " & Fmt.Int(n));
          Wr.PutText(stderr, "\n");
        END ReportAttemptToPlace;
        
      PROCEDURE ReportIdealPosition(n: Node; gp: IntPoint) =
        BEGIN
          Wr.PutText(stderr, "  ideal grid position of node " & Fmt.Int(n));
          Wr.PutText(stderr, " is ");
          PrintIntPoint(stderr, gp);
          Wr.PutText(stderr, "\n");
        END ReportIdealPosition;
        
      PROCEDURE ReportInsufficientHorRange(n: Node; horRange: RealRange) =
        BEGIN
          Wr.PutText(stderr, "  hor range for node ");
          Wr.PutText(stderr, Fmt.Int(n));
          Wr.PutText(stderr, " is ");
          PrintRealRange(stderr, horRange);
          Wr.PutText(stderr, " = ");
          WITH
            r = IntRange{
              lo := CEILING(horRange.lo/gStep),
              hi := FLOOR(horRange.hi/gStep)
            }
          DO
            PrintIntRange(stderr, r)
          END;
          Wr.PutText(stderr, " -- needs to unplace neighbors\n")
        END ReportInsufficientHorRange;
      
      PROCEDURE ReportFailureToPlace(n: Node) =
        BEGIN
          Wr.PutText(stderr, "  attempt to place node ");
          Wr.PutText(stderr, Fmt.Int(n));
          Wr.PutText(stderr, "  failed\n");
        END ReportFailureToPlace;
      
      PROCEDURE ReportPlacement(n: Node; gp: IntPoint) =
        BEGIN
          Wr.PutText(stderr, "  placing node " & Fmt.Int(n));
          Wr.PutText(stderr, " at ");
          PrintIntPoint(stderr, gp);
          Wr.PutText(stderr, "\n");
        END ReportPlacement;
        
      PROCEDURE ReportUnplacement(n: Node; gp: IntPoint) =
        BEGIN
          Wr.PutText(stderr, "  unplacing node " & Fmt.Int(n));
          Wr.PutText(stderr, " from ");
          PrintIntPoint(stderr, gp);
          Wr.PutText(stderr, "\n");
        END ReportUnplacement;

      PROCEDURE ComputeHorRange(n: Node): RealRange =
      (*
        Computes the min and max /x/ coordinate of node /n/,
        taking into account the current positions of its placed neighbors. 
        *)
        VAR xmin: REAL := 0.0;      (* Lower bound for /x/ *)
            xmax: REAL := gDims[0]; (* Upper bound for /x/ *)
            ri: REAL := nodeData[n].ri;
            ro: REAL := nodeData[n].ro;

        PROCEDURE ProcessSucc (* : Reduced.ArcAction *) (
            <*UNUSED*> i: NAT; 
            a: Arc;
          ) RAISES {} =
        (*
          Processes an outgoing edge of node /n/, updating /xmax/
          and accumulating the /x/ and /y/ sum. *)
          BEGIN
            WITH succ = node[a.dest] DO
              IF succ # NullNode AND placed[succ] THEN
                WITH succD = nodeData[succ] DO
                  <* ASSERT succD.gp[0] # Nowhere[0] *>
                  <* ASSERT succD.gp[1] # Nowhere[1] *>
                  xmax := MIN(xmax, XYFromGP(succD.gp[0]) - ro - succD.ri);
                END
              END
            END
          END ProcessSucc;

        PROCEDURE ProcessPred (* : Reduced.ArcAction *) (
            <*UNUSED*> i: NAT; 
            a: Arc;
          ) RAISES {} =
          BEGIN
            WITH pred = node[a.dest] DO
              IF pred # NullNode AND placed[pred] THEN
                WITH predD = nodeData[pred] DO
                  <* ASSERT predD.gp[0] # Nowhere[0] *>
                  <* ASSERT predD.gp[1] # Nowhere[1] *>
                  xmin := MAX(xmin, XYFromGP(predD.gp[0]) + ri + predD.ro);
                END
              END
            END
          END ProcessPred;

        BEGIN
          aut.EnumInArcs(nodeData[n].state, action := ProcessPred);
          aut.EnumOutArcs(nodeData[n].state, action := ProcessSucc);
          WITH 
            r = RealRange{xmin, xmax}
          DO
            <* ASSERT r.lo >= 0.0 *>
            <* ASSERT r.hi <= gDims[0] *>
            RETURN r
          END
        END ComputeHorRange;

      PROCEDURE UnplaceNearestNeighbors(n: Node; horRange: RealRange) =
      (*
        Given a node /n/ whose horizontal range /horRange/ is too
        narrow (or empty), attempts to unplace the placed predecessors 
        and successors of /n/ that determine that horizontal range, unless 
        they are fixed.
        Bombs out if the unplacements done are not enough to widen
        the horizontal range.
        *)
        VAR ri: REAL := nodeData[n].ri;
            ro: REAL := nodeData[n].ro;
            loPinned: BOOL := (horRange.lo <= 0.0);
            hiPinned: BOOL := (horRange.hi >= gDims[0]);

        PROCEDURE ProcessSucc (* : Reduced.ArcAction *) (
            <*UNUSED*> i: NAT; 
            a: Arc;
          ) RAISES {} =
        (*
          Processes an outgoing edge of node /n/, updating /xmax/
          and accumulating the /x/ and /y/ sum. *)
          BEGIN
            WITH succ = node[a.dest] DO
              IF succ # NullNode AND placed[succ] THEN
                WITH succD = nodeData[succ] DO
                  <* ASSERT succD.gp[0] # Nowhere[0] *>
                  <* ASSERT succD.gp[1] # Nowhere[1] *>
                  IF horRange.hi >= XYFromGP(succD.gp[0]) - ro - succD.ri THEN
                    IF fixed[succ] THEN
                      hiPinned := TRUE
                    ELSE
                      UnplaceNode(succ);
                      PutUnplacedNodeInQueue(succ)
                    END
                  END
                END
              END
            END
          END ProcessSucc;

        PROCEDURE ProcessPred (* : Reduced.ArcAction *) (
            <*UNUSED*> i: NAT; 
            a: Arc;
          ) RAISES {} =
          BEGIN
            WITH pred = node[a.dest] DO
              IF pred # NullNode AND placed[pred] THEN
                WITH predD = nodeData[pred] DO
                  <* ASSERT predD.gp[0] # Nowhere[0] *>
                  <* ASSERT predD.gp[1] # Nowhere[1] *>
                  IF horRange.lo <= XYFromGP(predD.gp[0]) + ri + predD.ro THEN
                    IF fixed[pred] THEN
                      loPinned := TRUE
                    ELSE
                      UnplaceNode(pred);
                      PutUnplacedNodeInQueue(pred)
                    END
                  END
                END
              END
            END
          END ProcessPred;

        BEGIN
          aut.EnumInArcs(nodeData[n].state, action := ProcessPred);
          aut.EnumOutArcs(nodeData[n].state, action := ProcessSucc);
          IF hiPinned AND loPinned THEN
            Wr.PutText(stderr, "** Not enough space between fixed nodes **\n");
            RAISE BadParameters
          END
        END UnplaceNearestNeighbors;

      PROCEDURE IdealPosition(n: Node): IntPoint =
      (*
        Computes the ideal grid indices of node /n/ from the 
        coordinates of its placed neighbors (predecessors and successors),
        or from the state number if it has no placed neighbors. *)

        VAR x: REAL := 0.0;
            y: REAL := 0.0;
            dydx: REAL := 0.0;
            nEdges: NAT := 0;

            ri: REAL := nodeData[n].ri;
            ro: REAL := nodeData[n].ro;

        PROCEDURE ProcessSucc (* : Reduced.ArcAction *) (i: NAT; a: Arc) RAISES {} =
        (*
          Processes an outgoing edge of node /n/, updating /xmax/
          and accumulating the /x/ and /y/ sum. *)
          BEGIN
            WITH succ = node[a.dest] DO
              IF succ # NullNode AND placed[succ] THEN
                WITH succD = nodeData[succ] DO
                  <* ASSERT succD.gp[0] # Nowhere[0] *>
                  <* ASSERT succD.gp[1] # Nowhere[1] *>
                  x := x + XYFromGP(succD.gp[0]) - ro - succD.ri;
                  WITH
                    slope = 0.5 * ComputeInitialArcSlope(n, i)
                  DO
                    y := y + XYFromGP(succD.gp[1]) - slope * XYFromGP(succD.gp[0]);
                    dydx := dydx + slope
                  END;
                  INC(nEdges);
                END
              END
            END
          END ProcessSucc;

        PROCEDURE ProcessPred (* : Reduced.ArcAction *) (i: NAT; a: Arc) RAISES {} =
          BEGIN
            WITH pred = node[a.dest] DO
              IF pred # NullNode AND placed[pred] THEN
                WITH predD = nodeData[pred] DO
                  <* ASSERT predD.gp[0] # Nowhere[0] *>
                  <* ASSERT predD.gp[1] # Nowhere[1] *>
                  x := x + XYFromGP(predD.gp[0]) + ri + predD.ro;
                  WITH
                    slope = 0.5 * ComputeInitialArcSlope(pred, i)
                  DO
                    y := y + XYFromGP(predD.gp[1]) - slope * XYFromGP(predD.gp[0]);
                    dydx := dydx + slope
                  END;
                  INC(nEdges)
                END
              END
            END
          END ProcessPred;

        BEGIN
          aut.EnumInArcs(nodeData[n].state, action := ProcessPred);
          aut.EnumOutArcs(nodeData[n].state, action := ProcessSucc);

          IF nEdges = 0 THEN 
            RETURN IdealPositionFromStateNumber(nodeData[n].state)
          END;

          x := x / FLOAT(nEdges);

          dydx := dydx / FLOAT(nEdges);
          y := y / FLOAT(nEdges) + dydx * x;

          IF o.verbose THEN 
            Wr.PutText(stderr, "  ideal coordinates of node ");
            Wr.PutText(stderr, Fmt.Int(n));
            Wr.PutText(stderr, " are ");
            PrintRealPoint(stderr, RealPoint{x, y});
            Wr.PutText(stderr, "\n")
          END;
          RETURN IntPoint{GPFromXY(x), GPFromXY(y)}
        END IdealPosition;
        
      PROCEDURE ComputeInitialArcSlope(n: Node; i: NAT): REAL =
      (*
        Computes the starting slope of the /i/th arc out of node /n/.
        *)
        BEGIN
          WITH
            s = nodeData[n].state,
            odeg = aut.OutDeg(s)
          DO
            RETURN RTan(FLOAT(i+1)/FLOAT(odeg+2))
          END
        END ComputeInitialArcSlope;
        
      PROCEDURE ComputeBestPlacement(
          n: Node; 
          gpIdeal: IntPoint; 
          horRange: RealRange;
        ): IntPoint =
      (*
        Computes the best placement for node /n/, whose ideal placement 
        would be at /gpIdeal/, taking into acount the position of already 
        placed nodes.

        The /horRange/ is the range of valid horizontal coordinates for /n/.
        *)
        VAR
          gRange: IntRange := IntRange{
            lo := CEILING(horRange.lo/gStep),
            hi := FLOOR(horRange.hi/gStep)
          };
          gpBest: IntPoint;
          penaltyBest: REAL;
        BEGIN
          <* ASSERT gRange.lo <= gRange.hi *>
          gpBest[0] := MAX(gRange.lo, MIN(gRange.hi, gpIdeal[0]));
          gpBest[1] := MAX(0, MIN(gridN[1], gpIdeal[1]));
          penaltyBest := CollisionPenalty(n, gpBest);
          IF penaltyBest <= 0.0 THEN RETURN gpBest END;

          CONST
            MaxTry = 10;  (* Number or alternative placements to try *)
          VAR 
            dir: ARRAY [0..1] OF [-1..+1];  (* Quadrant of /gpIdeal/ rel. center *)
            size: Counts;      (* Size of target region where to throw each dart *)
            base: IntPoint;    (* Corner of target region nearest to center *)
            gp: IntPoint;
          BEGIN
            IF gpIdeal[0] > gridN[0] DIV 2 THEN dir[0] := +1 ELSE dir[0] := -1 END;
            IF gpIdeal[1] > gridN[1] DIV 2 THEN dir[1] := +1 ELSE dir[1] := -1 END;
            WITH
              iri = CEILING(nodeData[n].ri),
              iro = CEILING(nodeData[n].ro),
              iht = MAX(iri, iro)
            DO
              size[0] := 1 + iri + iro + 1;
              size[1] := 1 + 2*iht + 1;
            END;
            base[0] := gpIdeal[0] - size[0]*dir[0];
            base[1] := gpIdeal[1] - size[1]*dir[1];
            FOR try := 1 TO MaxTry DO
              FOR c := 0 TO 1 DO
                gp[c] := base[c] + dir[c] * Random.Subrange(NIL, 0, size[c]);
                gp[c] := MAX(gRange.lo, MIN(gRange.hi, gp[c]));
                base[c] := base[c] + 2 * dir[c]
              END;
              WITH p = CollisionPenalty(n, gp) DO
                IF o.verbose THEN
                  Wr.PutText(stderr, "    trying node ");
                  Wr.PutText(stderr, Fmt.Int(n));
                  Wr.PutText(stderr, " at ");
                  PrintIntPoint(stderr, gp);
                  Wr.PutText(stderr, " penalty = ");
                  Wr.PutText(stderr, Fmt.Real(p, 3));
                  Wr.PutText(stderr, "\n")
                END;
                IF p < penaltyBest THEN
                  penaltyBest := p;
                  gpBest := gp
                END
              END;
              IF penaltyBest <= 0.0 THEN RETURN gpBest END;
            END;
          END;
          IF penaltyBest >= ExcessivePenalty THEN
            RETURN Nowhere
          ELSE
            RETURN gpBest
          END
        END ComputeBestPlacement;
        
      CONST
        ExcessivePenalty = 1.0e20;  (* Cannot place there, period *)

      PROCEDURE CollisionPenalty(n: Node; gp: IntPoint): REAL =
      (*
        Penalty (due to colisions with already placed nodes)
        if we were to place node /n/ at position /gp/. 
        *)
        VAR 
          penalty: REAL := 0.0;
        BEGIN
          WITH 
            iri = CEILING(nodeData[n].ri),
            iro = CEILING(nodeData[n].ro),
            iht = MAX(iri, iro)
          DO
            FOR i0 := MAX(0, gp[0] - iri) TO MIN(gridN[0], gp[0] + iro) DO
              FOR i1 := MAX(0, gp[1] - iht) TO MIN(gridN[1], gp[1] + iht) DO
                WITH w = who[i0, i1] + 0 DO
                  IF w # NullNode THEN
                    WITH p = SingleCollisionPenalty(w) DO
                      IF p = ExcessivePenalty THEN 
                        RETURN ExcessivePenalty 
                      ELSE
                        penalty := penalty + p
                      END;
                    END
                  END
                END
              END
            END
          END;
          RETURN penalty
        END CollisionPenalty;
        
      PROCEDURE SingleCollisionPenalty(w: Node): REAL =
      (*
        Penalty for colliding with a placed node /w/ in one grid cell.
        Returns ExcessivePenalty if /w/ is permanently fixed. 
        *)
        BEGIN
          <* ASSERT w # NullNode *>
          IF fixed[w] THEN
            RETURN ExcessivePenalty
          ELSE
            RETURN FLOAT(when[w])
          END
        END SingleCollisionPenalty;
        
      PROCEDURE UnplaceConflictingNodes(n: Node; gp: IntPoint) =
      (*
        Unplaces all placed nodes that collide with node /n/ 
        if the latter were to be placed at position /gp/,
        and puts them in the umplaced node queue: 
        *)
        BEGIN
          WITH 
            iri = CEILING(nodeData[n].ri),
            iro = CEILING(nodeData[n].ro),
            iht = MAX(iri, iro)
          DO
            FOR i0 := MAX(0, gp[0] - iri) TO MIN(gridN[0], gp[0] + iro) DO
              FOR i1 := MAX(0, gp[1] - iht) TO MIN(gridN[1], gp[1] + iht) DO
                WITH w = who[i0, i1] + 0 DO
                  IF w # NullNode THEN
                    UnplaceNode(w);
                    PutUnplacedNodeInQueue(w)
                  END
                END
              END
            END
          END
        END UnplaceConflictingNodes;
        
      PROCEDURE UnplaceNode(n: Node) =
      (*
        Unplaces the placed node /n/, and erases it from the /who/ map.
        *)
        BEGIN
          <* ASSERT placed[n] *>
          placed[n] := FALSE;
          WITH 
            gp = nodeData[n].gp,
            iri = CEILING(nodeData[n].ri),
            iro = CEILING(nodeData[n].ro),
            iht = MAX(iri, iro)
          DO
            (* Print unplacement report *)
            ReportUnplacement(n, gp);
            
            (* Now erase it from /who/ map: *)
            FOR i0 := MAX(0, gp[0] - iri) TO MIN(gridN[0], gp[0] + iro) DO
              FOR i1 := MAX(0, gp[1] - iht) TO MIN(gridN[1], gp[1] + iht) DO
                WITH w = who[i0, i1] DO
                  <* ASSERT w = NullNode OR w = n *>
                  w := NullNode
                END
              END
            END
          END
        END UnplaceNode;
        
      PROCEDURE PlaceNode(n: Node; gp: IntPoint; now: NAT) =
      (*
        Places the unplaced node /n/, and paints it on the /who/ map.
        *)
        BEGIN
          <* ASSERT NOT placed[n] *>
          placed[n] := TRUE;
          when[n] := now;
          WITH 
            iri = CEILING(nodeData[n].ri),
            iro = CEILING(nodeData[n].ro),
            iht = MAX(iri, iro)
          DO
            nodeData[n].gp := gp;
            FOR i0 := MAX(0, gp[0] - iri) TO MIN(gridN[0], gp[0] + iro) DO
              FOR i1 := MAX(0, gp[1] - iht) TO MIN(gridN[1], gp[1] + iht) DO
                WITH w = who[i0, i1] DO
                  <* ASSERT w = NullNode *>
                  w := n
                END
              END
            END
          END
        END PlaceNode;
        
      PROCEDURE FindNode(n: Node; READONLY a: ARRAY OF Node): NAT =
      (*
        Returns the position of /n/ in /a/, or NUMBER(a) if /n/ is not there. 
        *)
        BEGIN
          FOR i := 0 TO LAST(a) DO
            IF a[i] = n THEN RETURN i END
          END;
          RETURN NUMBER(a)
        END FindNode;

      PROCEDURE XYFromGP(i: INT): REAL =
        BEGIN
          RETURN FLOAT(i) * gStep
        END XYFromGP;

      PROCEDURE GPFromXY(READONLY x: REAL): INT =
        BEGIN
          RETURN ROUND(x / gStep)
        END GPFromXY;

      PROCEDURE PrintGridData() =
        BEGIN
          Wr.PutText(stderr, "Grid dimensions = ");
          Wr.PutText(stderr, Fmt.Real(gDims[0], 3, Fmt.Style.Flo));
          Wr.PutText(stderr, "m x ");
          Wr.PutText(stderr, Fmt.Real(gDims[1], 3, Fmt.Style.Flo));
          Wr.PutText(stderr, "m\n");
          Wr.PutText(stderr, "Grid step = ");
          Wr.PutText(stderr, Fmt.Real(gStep, 3, Fmt.Style.Flo));
          Wr.PutText(stderr, "m\n");
          Wr.PutText(stderr, "Grid size = ");
          Wr.PutText(stderr, Fmt.Int(gridN[0]));
          Wr.PutText(stderr, " x ");
          Wr.PutText(stderr, Fmt.Int(gridN[1]));
          Wr.PutText(stderr, "\n");
        END PrintGridData;

      BEGIN
        IF o.verbose THEN PrintGridData() END;
        PlaceNodes(effort := 3*maxNode);
      END;
      RETURN drawing
    END
  END DrawAutomaton;
  
EXCEPTION
  BadParameters;    (* Bad input parameters *)
  BadChar(CHAR);
  BadLetter(Basics.Letter);

PROCEDURE PrintLetter(wr: Wr.T; letter: Basics.Letter) =
  BEGIN
    TRY
      Wr.PutChar(wr, LetterToChar(letter))
    EXCEPT
    | BadLetter(bad) => 
        Wr.PutChar(wr, '(');
        WrNat(wr, ORD(bad));
        Wr.PutChar(wr, ')');
    END;
  END PrintLetter;

PROCEDURE CharToLetter(c: CHAR): Reduced.Letter RAISES {BadChar} =
  BEGIN
    IF c IN Char.Controls THEN RAISE BadChar(c) END;
    RETURN ORD(c)
  END CharToLetter;
  
PROCEDURE LetterToChar(letter: Basics.Letter): CHAR RAISES {BadLetter} =
  BEGIN
    IF letter = 0 THEN RAISE BadLetter(letter) END;
    WITH c = VAL(letter, CHAR) DO
      IF c IN Char.Controls THEN RAISE BadLetter(letter) END;
      RETURN c
    END
  END LetterToChar;

PROCEDURE GetOptions(): Options =
  VAR o: Options;
  BEGIN
    ParseParams.BeginParsing(stderr);
      
      ParseParams.GetKeyword("-load");
      o.loadFileName := ParseParams.GetNext();
      
      IF ParseParams.KeywordPresent("-print") THEN
        o.printFileName := ParseParams.GetNext()
      ELSE
        o.printFileName := ""
      END;
      
      IF ParseParams.KeywordPresent("-ps") THEN
        o.psFileName := ParseParams.GetNext()
      ELSE
        o.psFileName := ""
      END;
      
      o.verbose := ParseParams.KeywordPresent("-verbose");

      IF ParseParams.KeywordPresent("-gridDims") THEN
        o.gridDims[0] := ParseParams.GetNextReal(0.01, 100.0);
        o.gridDims[1] := ParseParams.GetNextReal(0.01, 100.0);
      ELSE
        o.gridDims[0] := 0.15;
        o.gridDims[1] := 0.12;
      END;
      
      IF ParseParams.KeywordPresent("-arcSpacing") THEN
        o.arcSpacing := ParseParams.GetNextReal(0.001, 0.010);
      ELSE
        o.arcSpacing := 0.003;
      END;
      
      WITH md = MAX(o.gridDims[0], o.gridDims[1]) DO      
        IF ParseParams.KeywordPresent("-gridStep") THEN
          o.gridStep := ParseParams.GetNextReal(md/FLOAT(MaxGridN), md/FLOAT(MinGridN));
        ELSE
          o.gridStep := md/FLOAT(MaxGridN)
        END;
      END;

    ParseParams.EndParsing();
    RETURN o
  END GetOptions;

PROCEDURE PrintIntPoint(wr: Wr.T; gp: IntPoint) =
  BEGIN
    Wr.PutChar(wr, '<');
    Wr.PutText(wr, Fmt.Int(gp[0]));
    Wr.PutChar(wr, ',');
    Wr.PutChar(wr, ' ');
    Wr.PutText(wr, Fmt.Int(gp[1]));
    Wr.PutChar(wr, '>');
  END PrintIntPoint;

PROCEDURE PrintRealPoint(wr: Wr.T; rp: RealPoint) =
  BEGIN
    Wr.PutChar(wr, '(');
    Wr.PutText(wr, Fmt.Real(rp[0], 3, Fmt.Style.Flo));
    Wr.PutChar(wr, ',');
    Wr.PutChar(wr, ' ');
    Wr.PutText(wr, Fmt.Real(rp[1], 3, Fmt.Style.Flo));
    Wr.PutChar(wr, ')');
  END PrintRealPoint;

PROCEDURE PrintIntRange(wr: Wr.T; r: IntRange) =
  BEGIN
    Wr.PutChar(wr, '[');
    Wr.PutText(wr, Fmt.Int(r.lo));
    Wr.PutChar(wr, '.');
    Wr.PutChar(wr, '.');
    Wr.PutText(wr, Fmt.Int(r.hi));
    Wr.PutChar(wr, ']');
  END PrintIntRange;

PROCEDURE PrintRealRange(wr: Wr.T; r: RealRange) =
  BEGIN
    Wr.PutChar(wr, '[');
    Wr.PutText(wr, Fmt.Real(r.lo, 3, Fmt.Style.Flo));
    Wr.PutChar(wr, '_');
    Wr.PutChar(wr, '_');
    Wr.PutText(wr, Fmt.Real(r.hi, 3, Fmt.Style.Flo));
    Wr.PutChar(wr, ']');
  END PrintRealRange;

PROCEDURE RTan(u: REAL): REAL =
(*
  Tangent of u * Pi/2 *)
  BEGIN
    RETURN FLOAT(Math.tan(LONGFLOAT(u*3.1415926)))
  END RTan;

BEGIN
  Main()
END Draw.
