(* Last edited on 1999-06-05 20:50:17 by stolfi *)

INTERFACE DAGF;

(* Implementation details and low-level operations for "DAG.T". *)
(* See the copyright and disclaimer note at the end of this file. *)

IMPORT DAG, BinGraph;
FROM DAG IMPORT Node;

TYPE
  T = DAG.T;

REVEAL
  DAG.Parent = BinGraph.T BRANDED OBJECT END;
    (*
      A "DAG.T" is actually a special case of "BinGraph.T",
      with some restrictions on the arcs.  The basic "DAG.T" methods
      are the basic "BinGraph.T" methods, with different names
      (and different mental pictures).
      
      REPRESENTATION CONVENTIONS: A "DAG.T" is actually a "BinGraph.T", with 
      the following correspondences:

          "Class(s) = LSym(s)"

          "HasArcs(s) = (RNode(s) # s)"  (i.e. "s" is a sink iff "RNode(s) = s")

          "Rest(s) = LNode(a)"

          "Last(s).dest = RNode(s)"

          "Last(s).label = RSym(s)"

      DAG INVARIANTS: A "DAG.T" also imposes the following invariants
      on the underlying "BinGraph.T":

        SINKREP:

           "(RNode(s) = s)  =>  (LNode(s) = s) AND (LSym(s) = 0)

        TOPORDER:

           "(RNode(s) # s)  =>  (LNode(s) < s) AND (RNode(s) < s)

        HOMOCLASS:

           "LSym(LNode(s)) = LSym(s)"
    *)

REVEAL
  DAG.T <: Private;

TYPE 
  Private = DAG.Public OBJECT
    METHODS

      Alloc(maxNodes: NAT): T;
        (*
          Expands the internal storage, if necessary, to accomodate up
          to "maxNodes" distinct nodes, including all nodes already in
          the "DAG.T".  May de-allocate storage if "NAlloc()"
          is much greater than "maxNodes".
          
          Does not change "nNodes", "epoch", or the contents 
          of any existing node.  Requires "maxNodes >= NNodes()". *)

      NAlloc(): NAT;
        (*
          Max number of nodes for which internal storage is currently allocated.
          Clients can create "NAlloc() - NNodes()" new nodes without raising "Full". *)

      AddSink(c: Symbol): Node RAISES {Full};
        (*
          Adds a new sink node to the graph, with given class. *)
          
      AddArc(rest: Node; label: Symbol; dest: Node): Node RAISES {Full};
        (*
          Adds a new node to the graph, with given "Rest" and "Last". *)
          
      Discard(s: Node): T;
        (*
          Removes from the "DAG.T" the node "s" and all higher-numbered
          nodes. *)
      
      (***************************************************************************)
      (* INPUT/OUTPUT                                                            *)
      (***************************************************************************)

      Load(rd: Rd.T): T;
        (*
          Discards all current nodes and arcs of the "DAG.T",
          and replaces them by data read from "rd" (which must be in
          the format used by "Dump"). *)

    END;
    
  NAT = CARDINAL;

END DAGF.

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
