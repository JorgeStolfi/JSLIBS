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

INTERFACE ADDiagram;

FROM ADTypes IMPORT 
  NAT, BOOL, IntPoint, IntBox, Counts, ArcControlPoints;

CONST
  MaxGridN = 2000;       (* Maximum number of positions in grid, each axis *)
  MinGridN = 20;         (* Minimum number of positions in grid, each axis *)
  MaxNodes = 256*256-1;  (* Maximum number of nodes that can be placed *)

TYPE
  T <: Public;

  Public = OBJECT

      nNodes: NAT;               (* Number of nodes, not counting NIL *)
      node: REF ARRAY OF Node;   (* Maps NodeLabel to Node. *)
      gridN: Counts;             (* Number of grid cells in X and Y. *)
      unplaced: NodeQueue;   (* Queue of unplaced nodes *)
      unhappy: NodeQueue;    (* Queue of placed but unhappy nodes *)

    METHODS
    
      (****************************************************************)
      (* NODES                                                        *)
      (****************************************************************)
      
      AddNode(
        label: NodeLabel; (* Node label *)
        marked: BOOL; (* TRUE to highlight the node *)
        ri: NAT;      (* Radius of input fan, in grid cells *)
        ro: NAT;      (* Radius of output fan, in grid cells *)
        ht: NAT;      (* Half-height, in grid cells *)
      ): Node;
      (* 
        Adds a new node to the diagram, returns it.
        The node is initially unplaced, and its range is 
        everywhere. *)

      EnumNodes(action: NodeAction) RAISES {Abort};
      (*
        Applies "action" to all nodes. *)
        
      EnumNodeCells(n: Node; action: NodeCellAction) RAISES {Abort};
      (*
        Applies "action" to all grid cells covered by the given node. *)
        
      (****************************************************************)
      (* ARCS                                                         *)
      (****************************************************************)

      AddArc(
        label: ArcLabel; (* Arc label *)
        org: Node;       (* Origin node *)
        dst: Node;       (* Destination node *)
      ): Arc;
      (* 
        Adds a new arc to the diagram, returns it. *)
        
      EnumArcs(action: ArcAction) RAISES {Abort};
      (*
        Applies "action" to all arcs. *)
        
      EnumArcCells(a: Arc; action: ArcCellAction) RAISES {Abort};
      (*
        Applies "action" to every grid cell that is entered by 
        (or is near to) the arc "a". (Note: the "action" must not itself
        call EnumArcCells on the same diagram!) *)
      
      (****************************************************************)
      (* NODE PLACEMENT                                               *)
      (****************************************************************)

      MakeNodeUnhappy(n: Node);
      (* 
        Marks node "n" as unhappy, and puts it in the unhappy queue 
        (if not already there). *)

      MoveNode(n: Node; gp: IntPoint);
      (*
        Changes the position of node "n" to "gp". 
        Requires that "n" be currently unplaced. *)
      
      UnplaceNode(n: Node);
      (*
        Unplaces the placed node "n", erases it and its arcs
        from the "who" map, and adds it to the unplaced node queue. *)

      PlaceNode(n: Node; gp: IntPoint; happy: BOOLEAN; now: NAT);
      (*
        Places the unplaced node "n" at "gp", and paints 
        it on the "who" map. Automatically unplaces all 
        placed nodes that overlap it. 
        If not "happy", puts it in the unhappy node queue. *)
      
    END;

TYPE
  Fragment = ROOT BRANDED OBJECT
       (* A fragment of a diagram: an Arc, Node, or a set thereof. *)
    METHODS 
      
      EnumElements (arcAction: ArcAction; nodeAction: NodeAction) RAISES {Abort};
      (*
        Enumerates all arcs and nodes in the fragment. *)
        
      AddElement (x: Element): Fragment;
      (*
        Returns a copy of "self" with the element "x" included.
        Returns "self" itself if "x" is already in "self". *)

      SubElement (x: Element): Fragment;
      (*
        Returns a copy of "self" with the element "x" excluded.
        Returns "self" itself if "x" is not in "self". *)
        
      Includes (x: Element): BOOL;
      (* 
        TRUE iff "self" includes element "x" *)
        
    END;
  
  Element <: Fragment;  (* An Arc or a Node *)

TYPE
  Node <: NodePublic;
  NodePublic = Element BRANDED OBJECT
      label: NodeLabel;   (* Node label *)
      placed: BOOL;       (* TRUE iff node is currently placed. *)
      happy: BOOL;        (* FALSE iff any neighbor has been placed since "when" *)
      marked: BOOL;       (* TRUE if node is to be highlighted in drawing *)
      gp: IntPoint;       (* Coordinates of node center (if currently placed). *)
      range: IntBox;      (* Allowed range for gp. *)
      when: NAT;          (* Generation when node was last placed (if placed). *)
      ht: NAT;            (* Nominal node height (in grid steps) is 2*ht + 1 *)
      tilt: REAL;         (* Correction to be applied to the ".dir" fields of arcs. *)
      s: ARRAY [0..1] OF NodeSideData; (* Data for the input (0) and output (1) side of node *)
    END;
    
  NodeLabel = [0..MaxNodes];   (* Node label (external node number). *)

  NodeSideData = RECORD   
      rad: NAT;      (* Nominal node radius (in grid steps). *)
      arcs: Arc;     (* First icoming/outgoing arc *)
    END;

TYPE
  Arc <: ArcPublic;
  ArcPublic = Element BRANDED OBJECT
      label: ArcLabel;    (* Arc label *)
      s: ARRAY [0..1] OF ArcEndData; (* Data for origin (0) and destination (1) nodes *)
    END;
    
  ArcLabel = CHAR;   (* Arc label (letter) *)
  
  ArcEndData = RECORD              
      node: Node; (* Node of origin/destination *)
      next: Arc;  (* Next arc with same origin/destination *)
      dir: REAL;  (* Raw departure/arrival angle, from -Pi/2 to +Pi/2 *)
    END;

VAR
  quiet:   BOOL;   (* If true, omits placed/unplaced messages. *)
  verbose: BOOL;   (* If true, causes various disgnostic info to be printed. *)

EXCEPTION Abort;

TYPE
  NodeAction = PROCEDURE (n: Node) RAISES {Abort};
  ArcAction =  PROCEDURE (a: Arc)  RAISES {Abort};

  NodeCellAction = PROCEDURE (n: Node; i0, i1: NAT; VAR cell: Fragment) RAISES {Abort};
  ArcCellAction =  PROCEDURE (a: Arc;  i0, i1: NAT; VAR cell: Fragment) RAISES {Abort};

PROCEDURE New(maxLabel: CARDINAL; maxNodes: CARDINAL; gridN: Counts): T;
  (* 
    Creates a new data record for the given number of nodes, 
    with a cell grid with the given number of cells in each dimension. *)
    
PROCEDURE NodeArea(n: Node): NAT;
  (*
    Area of node /n/, in grid cells (i.e. how many entries of /who/
    it would cover if placed). *)

PROCEDURE Degree(n: Node; side: [0..1]): NAT;
  (* 
    Indegree (side=0) or outdegree (side=1) of node "n". *)

PROCEDURE ComputeIdealPosition(t: T; n: Node; READONLY range: IntBox): IntPoint;
  (*
    Computes the best position of node "n" within the given "range",
    considering the gravitational energy of "gp" and the strain energies
    of all arcs incident to "n" (independently of whether their other ends 
    are placed or unplaced). Assumes arc directions at "n" are completely
    unconstrained. Does not consider collisions. *)

CONST
  ExcessivePenalty = 1.0e20;  (* Cannot place there, period *)

PROCEDURE PositionPenalty(t: T; n: Node; cutoff: REAL): REAL;
  (*
    Penalty for placing the unplaced node "n" at "n.gp",
    based only on arc strain (ignoring collisions).  
    
    Takes into account all neighbors, placed or unplaced.  Assumes the
    "tilt" of node "n", and the "dir" of every arc "a" incident to
    "n", have already been recomputed for its current position.
    
    May abort the computation, and return ExcessivePenalty, if it
    finds out that the true answer would be bigger than "cutoff". *)

PROCEDURE CollisionPenalty(t: T; n: Node; now: NAT; cutoff: REAL): REAL;
  (*
    Penalty, due to colisions with already placed nodes and arcs,
    if we were to place the unplaced node "n" at position "n.gp".
    
    Also includes the penalty for collisions of "n" with placed
    arcs, and collisions of the arcs that connect "n" to placed
    neighbors with already placed arcs and nodes.  Assumes the
    "tilt" of node "n", and the "dir" of every arc "a" incident to
    "n", have already been recomputed for its current position.
    
    May abort the computation, and return ExcessivePenalty, if it
    finds out that the true answer would be bigger than "cutoff". *)

PROCEDURE ComputeArcControlPoints(a: Arc): ArcControlPoints;
  (* 
    Computes the control points "cp" for drawing arc "a".  The arc 
    consists of a straight segment from cp[0] to cp[1], followed
    by a Bezier cubic arc with control points cp[1..4]. *)

(* NODE QUEUES *)

TYPE
  NodeQueue = REF RECORD
      node: REF ARRAY OF Node;
      front: NAT;   (* node[front] is first node of queue (if any; else front=rear). *)
      rear: NAT;    (* node[rear] is first empty slot in queue (there is always one). *)
    END;
    
PROCEDURE FirstNodeInQueue(queue: NodeQueue): Node;
  (*
    Returns the first node from the given queue,
    or NullNode if the queue is empty. *)
    
PROCEDURE CycleNodesInQueue(queue: NodeQueue);
  (*
    Pushes the node that is currently at the front of the queue to
    the rear of the queue. *)

PROCEDURE QueueLength(queue: NodeQueue): NAT;
  (* 
    Number of nodes still unplaced. *)
  
END ADDiagram.
