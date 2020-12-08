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

MODULE ADReport;

IMPORT Fmt, Text, Wr, Thread;
FROM Stdio IMPORT stderr;
FROM ADTypes IMPORT
  NAT, IntPoint, RealPoint, IntRange, IntBox, RealRange, Counts, Dims;
FROM ADDiagram IMPORT
  Arc, NodeLabel;

PROCEDURE ReportQueueStatus(nUnplaced, nUnhappy: CARDINAL) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(stderr, "  queues: " & Fmt.Int(nUnplaced) & " unplaced, ");
    Wr.PutText(stderr, Fmt.Int(nUnhappy) & " unhappy\n")
  END ReportQueueStatus;

PROCEDURE ReportSuccess() =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(stderr, "Success!\n")
  END ReportSuccess;

PROCEDURE ReportAttemptToPlace(lab: NodeLabel; placed: BOOLEAN) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutChar(stderr, '\n');
    Wr.PutText(stderr, "attempting to ");
    IF placed THEN
      Wr.PutText(stderr, "adjust placement of");
    ELSE
      Wr.PutText(stderr, "place");
    END;
    Wr.PutText(stderr, " node " & Fmt.Int(lab));
    Wr.PutText(stderr, "\n");
  END ReportAttemptToPlace;

PROCEDURE ReportIdealCoordinates(lab: NodeLabel; nPreds, nSuccs: NAT; p: RealPoint; ksum: REAL) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(stderr, "  ideal coordinates of node ");
    Wr.PutText(stderr, Fmt.Int(lab));
    Wr.PutText(stderr, " are ");
    PrintRealPoint(stderr, p);
    Wr.PutText(stderr, "\n");
    Wr.PutText(stderr, "    nPreds = " & Fmt.Int(nPreds) & "  nSuccs = " & Fmt.Int(nSuccs) & "\n");
    Wr.PutText(stderr, "    ksum = " & Fmt.Real(ksum) & "\n");
  END ReportIdealCoordinates;

PROCEDURE ReportIdealPosition(lab: NodeLabel; gp: IntPoint) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(stderr, "  ideal grid position of node " & Fmt.Int(lab));
    Wr.PutText(stderr, " is ");
    PrintIntPoint(stderr, gp);
    Wr.PutText(stderr, "\n");
  END ReportIdealPosition;

PROCEDURE ReportCurrentRange(lab: NodeLabel; range: IntBox) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(stderr, "  range for node ");
    Wr.PutText(stderr, Fmt.Int(lab));
    Wr.PutText(stderr, " is ");
    PrintIntRange(stderr, range[0]);
    Wr.PutText(stderr, " x ");
    PrintIntRange(stderr, range[1]);
    Wr.PutText(stderr, "\n")
  END ReportCurrentRange;

PROCEDURE ReportInsufficientHorRange(lab: NodeLabel; horRange: IntRange) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(stderr, "  hor range for node ");
    Wr.PutText(stderr, Fmt.Int(lab));
    Wr.PutText(stderr, " is ");
    PrintIntRange(stderr, horRange);
    Wr.PutText(stderr, " -- needs to unplace neighbors\n")
  END ReportInsufficientHorRange;

PROCEDURE ReportFailureToPlace(lab: NodeLabel) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(stderr, "  attempt to place node ");
    Wr.PutText(stderr, Fmt.Int(lab));
    Wr.PutText(stderr, "  failed\n");
  END ReportFailureToPlace;
  
PROCEDURE ReportTentativePlacement(
    lab: NodeLabel; 
    gp: IntPoint; 
    cp, dp: REAL; 
    kind: TEXT := "";
  ) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(stderr, "    node ");
    Wr.PutText(stderr, Fmt.Int(lab));
    Wr.PutText(stderr, " at ");
    PrintIntPoint(stderr, gp);
    Wr.PutText(stderr, " penalty = ");
    Wr.PutText(stderr, Fmt.Real(cp, 3));
    Wr.PutText(stderr, " + ");
    Wr.PutText(stderr, Fmt.Real(dp, 3));
    Wr.PutText(stderr, " = ");
    Wr.PutText(stderr, Fmt.Real(cp+dp, 3));
    IF NOT Text.Empty(kind) THEN 
      Wr.PutText(stderr, " (" & kind & ")");
    END;
    Wr.PutText(stderr, "\n")
  END ReportTentativePlacement;

PROCEDURE ReportEnergyComputation(
    lab: NodeLabel; 
    gp: IntPoint; 
    dirFixed: BOOLEAN;
    energy, xForce, yForce: REAL;
  ) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(stderr, "      node ");
    Wr.PutText(stderr, Fmt.Int(lab));
    Wr.PutText(stderr, " at ");
    PrintIntPoint(stderr, gp);
    Wr.PutText(stderr, " dirs = ");
    IF dirFixed THEN 
      Wr.PutText(stderr, " fix")
    ELSE
      Wr.PutText(stderr, " var")
    END;
    Wr.PutText(stderr, " energy = ");
    Wr.PutText(stderr, Fmt.Real(energy, 4));
    Wr.PutText(stderr, " force = ( ");
    Wr.PutText(stderr, Fmt.Real(xForce, 4));
    Wr.PutText(stderr, "  ");
    Wr.PutText(stderr, Fmt.Real(yForce, 4));
    Wr.PutText(stderr, " )\n")
  END ReportEnergyComputation;  

PROCEDURE ReportArcStrain(
    a: Arc;
    rx, ry, angm, angn: REAL;
    SE, dSE_dx, dSE_dy: REAL;
    TE, dTE_dx, dTE_dy: REAL;
    BE, dBE_dx, dBE_dy: REAL;
  ) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(stderr, "      arc ");
    PrintArc(stderr, a);
    Wr.PutText(stderr, "  r = ( ");
    Wr.PutText(stderr, Fmt.Real(rx, 4));
    Wr.PutText(stderr, "  ");
    Wr.PutText(stderr, Fmt.Real(ry, 4));
    Wr.PutText(stderr, " )");
    Wr.PutText(stderr, "  angm = ");
    Wr.PutText(stderr, Fmt.Real(angm, 4));
    Wr.PutText(stderr, "  angn = ");
    Wr.PutText(stderr, Fmt.Real(angn, 4));
    Wr.PutText(stderr, "\n");
    Wr.PutText(stderr, "        ");
    Wr.PutText(stderr, "  SE = ");
    Wr.PutText(stderr, Fmt.Real(SE, 4));
    Wr.PutText(stderr, " dSE/dp = ( ");
    Wr.PutText(stderr, Fmt.Real(dSE_dx, 4));
    Wr.PutText(stderr, "  ");
    Wr.PutText(stderr, Fmt.Real(dSE_dy, 4));
    Wr.PutText(stderr, " )");
    Wr.PutText(stderr, "\n");
    Wr.PutText(stderr, "        ");
    Wr.PutText(stderr, "  TE = ");
    Wr.PutText(stderr, Fmt.Real(TE, 4));
    Wr.PutText(stderr, " dTE/dp = ( ");
    Wr.PutText(stderr, Fmt.Real(dTE_dx, 4));
    Wr.PutText(stderr, "  ");
    Wr.PutText(stderr, Fmt.Real(dTE_dy, 4));
    Wr.PutText(stderr, " )");
    Wr.PutText(stderr, "\n");
    Wr.PutText(stderr, "        ");
    Wr.PutText(stderr, "  BE = ");
    Wr.PutText(stderr, Fmt.Real(BE, 4));
    Wr.PutText(stderr, " dBE/dp = ( ");
    Wr.PutText(stderr, Fmt.Real(dBE_dx, 4));
    Wr.PutText(stderr, "  ");
    Wr.PutText(stderr, Fmt.Real(dBE_dy, 4));
    Wr.PutText(stderr, " )");
    Wr.PutText(stderr, "\n");
  END ReportArcStrain;

PROCEDURE ReportArcArcCollision(a, b: Arc; i0, i1: INTEGER) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(stderr, "    arc-arc collision at ");
    PrintIntPoint(stderr, IntPoint{i0, i1});
    Wr.PutText(stderr, " a = ");
    PrintArc(stderr, a);
    Wr.PutText(stderr, " b = ");
    PrintArc(stderr, b);
    Wr.PutText(stderr, "\n");
  END ReportArcArcCollision;
  
PROCEDURE PrintArc(wr: Wr.T; a: Arc) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
   Wr.PutText(wr, "(");
   Wr.PutText(wr, Fmt.Int(a.s[0].node.label));
   Wr.PutText(wr, ", ");
   Wr.PutChar(wr, a.label);
   Wr.PutText(wr, ", ");
   Wr.PutText(wr, Fmt.Int(a.s[1].node.label));
   Wr.PutText(wr, ")");
  END PrintArc;

PROCEDURE ReportPlacement(lab: NodeLabel; gp: IntPoint) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(stderr, "  placing node " & Fmt.Int(lab));
    Wr.PutText(stderr, " at ");
    PrintIntPoint(stderr, gp);
    Wr.PutText(stderr, "\n");
  END ReportPlacement;

PROCEDURE ReportUnplacement(lab: NodeLabel; gp: IntPoint) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(stderr, "  unplacing node " & Fmt.Int(lab));
    Wr.PutText(stderr, " from ");
    PrintIntPoint(stderr, gp);
    Wr.PutText(stderr, "\n");
  END ReportUnplacement;

PROCEDURE ReportUnhappiness(lab: NodeLabel) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(stderr, "  node " & Fmt.Int(lab));
    Wr.PutText(stderr, " is now unhappy ");
    Wr.PutText(stderr, "\n");
  END ReportUnhappiness;

EXCEPTION FatalError;

PROCEDURE ErrorNotEnoughGridCells(totNodeCells: NAT) =
  <* FATAL Wr.Failure, Thread.Alerted, FatalError *>
  BEGIN
    Wr.PutText(stderr, "** Total node area = ");
    Wr.PutText(stderr, Fmt.Int(totNodeCells));
    Wr.PutText(stderr, " cells, too big for given grid **\n");
    RAISE FatalError
  END ErrorNotEnoughGridCells;

PROCEDURE ErrorNotEnoughSpace(lab: NodeLabel) =
  <* FATAL Wr.Failure, Thread.Alerted, FatalError *>
  BEGIN
    Wr.PutText(stderr, "** Not enough horizontal or vertical space for node ");
    Wr.PutText(stderr, Fmt.Int(lab));
    Wr.PutText(stderr, " **\n");
    RAISE FatalError
  END ErrorNotEnoughSpace;

PROCEDURE PrintIntPoint(wr: Wr.T; gp: IntPoint) =
   <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutChar(wr, '<');
    Wr.PutText(wr, Fmt.Int(gp[0]));
    Wr.PutChar(wr, ',');
    Wr.PutChar(wr, ' ');
    Wr.PutText(wr, Fmt.Int(gp[1]));
    Wr.PutChar(wr, '>');
  END PrintIntPoint;

PROCEDURE PrintRealPoint(wr: Wr.T; rp: RealPoint) =
   <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutChar(wr, '(');
    Wr.PutText(wr, Fmt.Real(rp[0], 3, Fmt.Style.Flo));
    Wr.PutChar(wr, ',');
    Wr.PutChar(wr, ' ');
    Wr.PutText(wr, Fmt.Real(rp[1], 3, Fmt.Style.Flo));
    Wr.PutChar(wr, ')');
  END PrintRealPoint;

PROCEDURE PrintIntRange(wr: Wr.T; r: IntRange) =
   <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutChar(wr, '[');
    Wr.PutText(wr, Fmt.Int(r.lo));
    Wr.PutChar(wr, '.');
    Wr.PutChar(wr, '.');
    Wr.PutText(wr, Fmt.Int(r.hi));
    Wr.PutChar(wr, ']');
  END PrintIntRange;

PROCEDURE PrintRealRange(wr: Wr.T; r: RealRange) =
   <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutChar(wr, '[');
    Wr.PutText(wr, Fmt.Real(r.lo, 3, Fmt.Style.Flo));
    Wr.PutChar(wr, '_');
    Wr.PutChar(wr, '_');
    Wr.PutText(wr, Fmt.Real(r.hi, 3, Fmt.Style.Flo));
    Wr.PutChar(wr, ']');
  END PrintRealRange;

PROCEDURE PrintGridData(gDims: Dims; gStep: REAL; gridN: Counts) =
  <* FATAL Wr.Failure, Thread.Alerted *>
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
END ADReport.
