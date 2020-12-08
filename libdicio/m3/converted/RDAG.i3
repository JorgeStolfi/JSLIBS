INTERFACE RDAG;

(* The RDAG data structure for reduced acyclic ordered labeled directed graphs. *)
(* See the copyright and disclaimer note at the end of this file. *)

(*
  An "RDAG.T" is a special case of a "DAG.T" (acyclic ordered labeled 
  directed graph) with non-redundant nodes.
  
  NODES, CLASSES, ARCS, LABELS: An "RDAG.T" inherits from a "DAG.T"
  the concepts of {\em node} (vertex) and {\em arc} (edge).  Each node
  has a {\em class}, and each edge has a {\em label} (both being
  "Symbol" values).  There may be more than one edge with same origin,
  destination, and label.
  
  The arcs out of a node "s" have a specific order; if "s" is not a
  sink, "Last(s)" gives the last arc in this order, and "Rest(s)"
  gives a node whose class is that of "s", and whose arc list is that
  of "s" minus "Last(s)".  See the "DAG" interface for more details.

  NODE UNIQUENESS: the implementation of an "RDAG.T" ensures also that
  nodes are uniquely represented, in the sense that no two distinct
  nodes have the same class and the same list of outgoing arcs
  (including multiplicity, order, labels, and destinations). 
  That is, the implementation ensures the following invariant:
  
  |   p # q  =>  Class(p) # Class(q) OR
  |              IsSink(p) # IsSink(q) OR
  |              NOT IsSink(p) AND ( Last(p) # Last(q) OR Rest(p) # Rest(q) )

  In particular, for each "Symbol" value "c", there is at most one
  sink node whose class is "c". 

*)

TYPE
  T <: Public;
  Public = DAG.T OBJECT
    METHODS

      AddArc(last: Arc; rest: Node): Node RAISES {Full};
      (*
        The unique node "s" such that "Last(s) = arc" and "Rest(s) =
        rest". (Creates the node if necessary).  Amortized cost:
        $O(1)$ time, $O(1)$ space.

        "AddArc" raises "Full" iff the currently allocated
        space for the BinGraph is exhausted (that is, if 
        "MaxNode() = MaxAllocNode()" and the node is not yet present).
        In this case the client can catch the exception
        and, for instance, invoke "Expand" or "Crunch" before 
        repeating the "AddArc". *)

      (***************************************************************************)
      (* STORAGE MANAGEMENT                                                      *)
      (***************************************************************************)

      MaxNode(): Node;
      (*
        Known nodes of this BinGraph are
        numbered from 0 to "MaxNode()". *)

      MaxAllocNode(): NAT;
      (*
        Implementation data: max node for which storage has been allocated.
        Clients can do at least "MaxAllocNode() - MaxNode()" calls to
        "AddArc" without raising "Full". *)

      Expand(newSize: NAT) RAISES {Full};
      (*
        Expands the internal storage to accomodate
        up to "newSize" distinct nodes, not counting NullNode
        but counting all nodes already in the BinGraph.
        Preserves all current nodes and their numbers. *)


      Crunch(VAR (*IO*) root: ARRAY OF Node);
      (*
        Discards all nodes not reachable from the given "root" nodes,
        squeezes the reachable ones together (preserving their order), and
        updates the "root" vector to reflect the new Node numbering.

        WARNING: By definition, "Crunch" generally changes the numbers
        of all nodes, and deletes some of them. The client must make
        sure that all important Node variables are part of the "root" vector.
        *)

    END;

TYPE
  ArcAction = PROCEDURE(len: NAT; org: Node; i: NAT; arc: Arc) RAISES {Abort, Skip};
  NodeAction = PROCEDURE(len: NAT; s: Node) RAISES {Abort, Skip};

(***************************************************************************)
(* CREATION                                                                *)
(***************************************************************************)

PROCEDURE New(size: NAT): T;
(*
  Creates a new BinGraph, initally with no arcs and a single node (NullNode).

  The "size" is the maximum number of distinct nodes (not counting
  "NullNode") that can be created by "AddArc" without raising "Full". *)

(***************************************************************************)
(* COPYING                                                                 *)
(***************************************************************************)

PROCEDURE Copy(
    from: T;
    to: T;
    s: Node;
    VAR (*IO*) map: REF ARRAY OF Node;
  ): Node RAISES {Full};
(*
  Copies into the "to" BinGraph all nodes of the "from" BinGraph
  that are reachable from "root".

  If the "map" argument is not NIL , "Copy" assumes that any node
  "t" with "map[t] # NullNode" has already been copied
  into the "to" BinGraph, and its number there is "map[t]".
  If "map" is NIL, "Copy" will allocate and return one of the
  appropriate size. In any case, "Copy" will set "map[t]"
  appropriately for every node that it copies, extending
  the array if necessary.

  "Copy" will maintain the uniqueneness invariant: when copying a node,
  it will create a new node in "to" only if there is no equivalent node
  there. (This will be true even if the existence of that node is not
  recorded in the client-given "map".)
  *)

(***************************************************************************)
(* IO                                                                      *)
(***************************************************************************)

PROCEDURE Dump(wr: Wr.T; BinGraph: T);
PROCEDURE Load(rd: Rd.T; minSize: NAT := 0): T;
(*
  "Dump" writes "BinGraph" to the given writer, in a  format that can be
  read back with "Load".  The output is in ASCII but not meant to be
  human-readable.

  "Dump" only writes the existing nodes, from "1" to "BinGraph.MaxNode()".
  The "BinGraph.T" returned by "Load" is just large enough to contain
  those nodes, or "minSize" nodes, whichever is larger. *)


  "Copy" will maintain the uniqueneness invariant: when copying a node,
  it will create a new node in "to" only if there is no equivalent node
  there. (This will be true even if the existence of that node is not
  recorded in the client-given "map".)
  *)

END RDAG.

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
