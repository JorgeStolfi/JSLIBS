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

INTERFACE ADReport;

IMPORT Wr;

FROM ADTypes IMPORT 
  NAT, Dims, Counts, RealRange, IntRange, IntBox, RealPoint, IntPoint;
FROM ADDiagram IMPORT
  Arc, NodeLabel;

(* Progress reports: *)

PROCEDURE ReportQueueStatus(nUnplaced, nUnhappy: CARDINAL);

PROCEDURE ReportSuccess();

PROCEDURE ReportAttemptToPlace(lab: NodeLabel; placed: BOOLEAN);

PROCEDURE ReportIdealCoordinates(lab: NodeLabel; nPreds, nSuccs: NAT; p: RealPoint; ksum: REAL);

PROCEDURE ReportIdealPosition(lab: NodeLabel; gp: IntPoint);

PROCEDURE ReportCurrentRange(lab: NodeLabel; range: IntBox);

PROCEDURE ReportInsufficientHorRange(lab: NodeLabel; horRange: IntRange);

PROCEDURE ReportFailureToPlace(lab: NodeLabel);

PROCEDURE ReportTentativePlacement(
    lab: NodeLabel; 
    gp: IntPoint; 
    cp, dp: REAL; 
    kind: TEXT := "";
  );

PROCEDURE ReportEnergyComputation(
    lab: NodeLabel; 
    gp: IntPoint; 
    dirFixed: BOOLEAN;
    energy, xForce, yForce: REAL;
  );

PROCEDURE ReportArcStrain(
    a: Arc;
    rx, ry, angm, angn: REAL;
    SE, dSE_dx, dSE_dy: REAL;
    TE, dTE_dx, dTE_dy: REAL;
    BE, dBE_dx, dBE_dy: REAL;
  );

PROCEDURE ReportArcArcCollision(a, b: Arc; i0, i1: INTEGER);

PROCEDURE ReportPlacement(lab: NodeLabel; gp: IntPoint);

PROCEDURE ReportUnplacement(lab: NodeLabel; gp: IntPoint);

PROCEDURE ReportUnhappiness(lab: NodeLabel);

(* Errors: *)

PROCEDURE ErrorNotEnoughSpace(lab: NodeLabel);

PROCEDURE ErrorNotEnoughGridCells(totNodeCells: NAT);

(* Low-level formatting: *)

PROCEDURE PrintIntPoint(wr: Wr.T; gp: IntPoint);

PROCEDURE PrintRealPoint(wr: Wr.T; rp: RealPoint);

PROCEDURE PrintIntRange(wr: Wr.T; r: IntRange);

PROCEDURE PrintRealRange(wr: Wr.T; r: RealRange);

PROCEDURE PrintGridData(gDims: Dims; gStep: REAL; gridN: Counts);

END ADReport.
