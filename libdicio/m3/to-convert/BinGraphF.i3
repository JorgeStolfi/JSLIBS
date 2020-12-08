(* Last edited on 1999-06-05 20:49:25 by stolfi *)

INTERFACE BinGraphF;

(* Implementation details and low-level operations for "BinGraph.T". *)
(* See the copyright and disclaimer note at the end of this file. *)

IMPORT BinGraph, Word, Wr, Rd;
FROM Basics IMPORT Full, NAT;

CONST
  NodePBits = 24;    (* Num bits to use for packed node numbers *)
  NodePMax = Word.Shift(1, NodePBits) - 1; (* Maximum packed node number *)

TYPE
  T = BinGraph.T;
  Node = BinGraph.Node;
  Symbol = BinGraph.Symbol;

  NodeP = BITS NodePBits FOR [0..NodePMax];

  Entry = RECORD  (* Data for one node/link *)
      lSym: Symbol;      (* Label of "L" link. *)
      lNode: NodeP;      (* Destination of "L" link. *)
      rSym: Symbol;      (* Label of "R" link. *)
      rNode: NodeP;      (* Destination of "R" link. *)
    END;

  Entries = ARRAY OF Entry;

REVEAL
  T <: Private;

TYPE 
  Private = BinGraph.Public OBJECT
      nNodes: NAT := 0;
      e: REF Entries;
        (*
          The existing nodes are "e[0..nNodes-1]". *)
          
      epoch: NAT := 0;
        (*
          A counter that can be used TO detect changes to the graph. 
          It is incremented every time an existing node is
          removed or modified (but not when new nodes are added).  
          Thus, if "s" is a valid node number, it will remain valid,
          and the graph rooted at "s" will remain unchanged, as long
          as "epoch" is unchanged.  Moreover, if "epoch" and "nNodes"
          are unchanged, then the whole graph is unchanged. *)

    METHODS
      
      Alloc(maxNodes: NAT): T;
        (*
          Expands the internal storage, if necessary, to accomodate up
          to "maxNodes" distinct nodes, including all nodes already in
          the "BinGraph.T".  May de-allocate storage if "NAlloc()"
          is much greater than "maxNodes".
          
          Does not change "NNodes", "Epoch", or the contents 
          of any existing node.  Requires "maxNodes >= NNodes()". *)

      NAlloc(): NAT;
        (*
          Max number of nodes for which internal storage is currently allocated.
          Clients can create "NAlloc() - NNodes()" new nodes without raising "Full". *)

      Set(s: Node; lSym, rSym: Symbol; lNode, rNode: Node);
        (*
          Modifies the data of node "s".  Requires that "s", "lNode",
          and "rNode" be less than "NNodes()".  Increments the "epoch"
          counter (if the data actually changed).  
          Cost: $O(1)$ time, 0 space. *)
      
      Add(lSym, rSym: Symbol; lNode, rNode: Node): Node RAISES {Full};
        (*
          Creates a new node with given labels and outgoing links.
          The new node's number (returned as a result) is equal to the
          value of "NNodes()" before the call.  Requires "lSym <=
          NNodes()" and "rSym <= NNodes()".  Increments "nNodes"
          but leaves "epoch" unchanged. 
          Cost: $O(1)$ time, $O(1)$ space.

          "Make" raises "Full" iff the currently allocated space for
          the "BinGraph.T" is exhausted (that is, "NNodes() = NAlloc()").
          In this case the client may catch the exception and, for instance,
          invoke "Expand" or "Crunch" before repeating the "Set". *)

      Clear(): T;
        (*
          Discards all nodes. Increments the "epoch" counter 
          (if there were any nodes). *)
      
      Copy(from: T; s: Node; VAR (*IO*) map: ARRAY OF Node): Node RAISES {Full};
        (*
          Copies into self every node "t" of the "from" graph that is
          reachable from "s".  
          
          If "map[t] # NoNode" "Copy" assumes that "t" has already been
          copied to node "map[t]" of self; otherwise "Copy" finds/creates
          an isomorphic node "t'" in self, and sets "map[t] := t'".
          May increment "nNodes" but does not modify the "epoch". *)
      
      Crunch(READONLY root: ARRAY OF Node; VAR (*OUT*) map: ARRAY OF Node): T;
        (*
          Discards all nodes not reachable from the given "root" nodes,
          squeezes the reachable ones together (preserving their order), 
          and returns the mapping between new and old node number in "map".
          Specifically, if "map[s] = NoNode" then node "s" was deleted,
          else its new number is "map[s]".  Increments the "epoch"
          (if any node has moved).

          Don't forget to update the "root" nodes as specified by "map". *)
          
      (***************************************************************************)
      (* INPUT/OUTPUT                                                            *)
      (***************************************************************************)

      Load(rd: Rd.T): T;
        (*
          Discards any existing nodes, links, and comment test of the graph,
          and replaces them by data read from "rd".
          
          The format of "rd" must be either that of "Dump",
          or the old "DAG.Dump" format of 96-11-16. *)

    END;
    
(***************************************************************************)
(* LOW_LEVEL INPUT/OUTPUT TOOLS                                            *)
(***************************************************************************)

PROCEDURE DumpBody(g: T; wr: Wr.T);
  (*
    Prints the nodes and links in the same format as "Dump", but 
    omitting the file header, footer, and comment text. *)

PROCEDURE LoadBody(g: T; rd: Rd.T);
  (*
    Discards all current nodes and links of "g", and replaces them by
    data read from "rd" (which must be in the format used by "Dump",
    minus the file header, comments, and footer). Preserves the comment
    text of "g", and tries to reuse storage if possible. *)

PROCEDURE LoadOldDAGBody(g: T; rd: T);
  (*
    Compatibility tool: similar to "LoadBody", but expects data in the format 
    created by the old version of "DAG.Dump" (format of "91-11-16),
    minus the file header, comments, and footer.  Automatically creates the
    "0" node (which was omitted in the old "DAG.Dump" format). Preserves the comment
    text of "g", and tries to reuse storage if possible. *)

END BinGraphF.

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
