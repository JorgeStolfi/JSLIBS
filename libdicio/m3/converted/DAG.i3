(* Last edited on 1999-06-05 20:49:41 by stolfi *)

INTERFACE DAG;

(* The DAG data structure for acyclic ordered labeled directed graphs. *)
(* See the copyright and disclaimer note at the end of this file. *)

(*
  A DAG.T is a finite directed acyclic graph, with labeled vertices
  and ordered labeled edges. Parallel edges are allowed.  
  
  NODES AND ARCS: The vertices and edges of the graph are called 
  {\em nodes} and {\em arcs} herein.
  
  CLASSES: Every node "s" of a "DAG.T" is labeled with a "Symbol",
  denoted by "Class(s)".
  
  ARC LABELS: Every arc has an {\em origin} node, a {\em destination}
  node, and a {\em label} (a "Symbol").

  SINK NODES: If the graph has any nodes at all, it must have at least 
  one ``sink'' (a node with no outgoing arcs).  The method call
  "IsSink(s)" returns "TRUE" iff node "s" is a sink.
  
  INTRINSIC ARC ORDER: for every node "s" that is not a sink, there
  is one distinguished arc "Last(s)" out of "s", and another node
  "Rest(s)", such that the set of arcs leaving the latter are those
  leaving "s", minus the arc "Last(s)".
  
  The methods "Last" and "Rest" define the ``intrinsic ordering'' for the
  arcs out "s".  Note that "Rest" moves {\em backwards} in this ordering.
    
  SUBNODES: by definition, a "subnode" of a node "s" of the DAG is a node
  "t" whose outgoing arcs are an initial prefix of the list of arcs
  out of "s", in intrinsic order. In other words, "t" is a subnode of 
  "s" if it can be obtained from "s" by zero or more applications of the
  "Rest" method.
  
  All subnodes of a node have the same class as that node, i.e.
  "Class(Rest(s)) = Class(s)" whenever "Rest(s)" is defined.  In
  particular, every "s" node has a subnode that is a sink of the
  same class as "s".

  TOPOLOGICAL ORDERING: The nodes in a "DAG.T" are
  identified by numbers, assigned as they are created.  The numerical
  order of the nodes is compatible with the DAG arcs and the
  subnode relation; that is, every node is greater (numerically)
  than any of its proper subnodes, and of the destination of any of
  its outgoing arcs.
  
  It follows that node 0, if it exists, is always a sink node.

  STANDARD PATH ORDERS: The intrinsic ordering of the arcs out of each node
  induces a ``standard'' ordering of the set of all paths out of a fixed node "s".
  In this ordering, all paths that leave "s" through the arc "Last(s)" occur after 
  those that leave through other arcs. 
  
  There are actually three standard path orderings, differing on where
  the paths that {\em end} at "s" appear: either before all the paths
  that {\em go through} "s" (``pre-order''), or after them
  (``end-order''), or both before and after them (``double order'').
  
*)

IMPORT Wr;
IMPORT Basics, Commented, BinGraph;
FROM Basics IMPORT NAT, BOOL;
FROM Basics IMPORT Skip, Abort;

TYPE
  Symbol = Basics.Symbol;
  String = Basics.String;
  Node = NAT;
  
  ArcData = RECORD
      label: Symbol;
      dest: Node
    END;
    
CONST
  NoNode = BinGraph.NoNode;

TYPE
  T <: Public;

TYPE
  Parent <: Commented.T;
  Public = Parent OBJECT
    METHODS

      NNodes(): NAT;
        (*
          Number of nodes in the automaton. *)

      Class(s: Node): Symbol;
        (*
          The node's class. *)
      
      IsSink(s: Node): BOOL;
        (*
          TRUE if "s" has one or more outgoing arcs; FALSE if "s" is a sink.
          Cost: $O(1)" time, 0 space. *)
      
      Last(s: Node): ArcData;
        (*
          The last arc out of node "s". Requires "IsSink(s)".
          Cost: $O(1)$ time, 0 space. *)

      Rest(s: Node): Node;
        (*
          A node "t" such that "Class(t) = Class(s)", and whose
          outgoing arcs are those of "s" minus "Last(s)".
          Requires "IsSink(s)".  Cost: $O(1)$ time, 0 space. *)

      (******************************************************************)
      (* DERIVED METHODS                                                *)
      (******************************************************************)

      (*
        These methods could be ordinary procedures, since
        they can be defined in terms of the primitives above;
        they are declared here as methods to reduce client confusion. (?)
      *)

      OutDeg(s: Node): NAT;
        (*
          Number of arcs out of node "s" (0 iff "s" is a sink).
          Cost: $O("OutDeg(s)")$ time, 0 space. *)

      First(s: Node): ArcData;
        (*
          The first arc out of node "s". Requires "IsSink(s)".
          Cost: $O("OutDeg(s)")$ time, 0 space. *)

      (*******************************************************************)
      (* ENUMERATION                                                     *)
      (*******************************************************************)

      EnumPaths(
          s: Node;
          enter: NodeAction := NIL;
          push, pop: ArcAction := NIL;
          exit: NodeAction := NIL;
        ) RAISES {Abort};
        (*
          Enumerates a set of paths in the DAG, starting from node "s",
          in depth-first order.  The set of paths is defined by the
          client-provided procedures "enter", "push", "pop", and "exit".

          The enumeration algorithm can be described in terms of a conceptual
          "current path" that initially contains just the node "s",
          and grows or shrinks one arc at a time.  "EnumPaths" calls

              "enter(len,o,c)" whenever the current path has length "len"
                  and has just reached node "o", whose class is "c".

              "push(len,o,c,i,r,d)" whenever the current path has "len" arcs,
                  ends at node "o", whose class is "c", and is being extended
                  with an arc labeled "r" to node "d", which is the "i"th arc
                  out of "o" (counting from "0 = First(o)").

              "pop(len,o,c,i,r,d)" whenever the current path has "len+1" arcs
                  and EnumPaths is about to remove its last arc, with label "r"
                  and destination "d", which is the "ith" arc out of "o",
                  whose class is "c".

              "exit(len,o,c)" whenever the current path has "len" arcs and
                  ends in node "o", whose class is "c", and no extensions
                  of it remain to be enumerated.

          Any of the four client actions can stop the enumeration by raising "Abort".

          Unless the enumeration is aborted, every call to "enter" will be followed
          eventually by a matching call to "exit", and every call to "push"
          will be followed eventually by a matching call to "pop".  The typical
          action pattern for a generic node "o" with "n" outgoing arcs is

          |     enter(len,o,c)
          |       push(len,o,c,0,x,d)
          |         ...
          |       pop(len,o,c,0,x,d)
          |       push(len,o,c,1,y,e)
          |         ...
          |       pop(len,o,c,1,y,e)
          |       ...
          |       push(len,o,c,n-1,z,f)
          |         ...
          |       pop(len,o,c,n-1,z,f)
          |     exit(len,o,c)

          In particular, "EnumPaths" will call "enter(0,s,Class(s))" at the very
          beginning, and (if not aborted) will call "exit(0,s,Class(s))" at the very end.

          Note that the same node "o" may be "enter"ed and "exit"ed many, many times.

          The client can prune branches of the path tree by raising the
          "Skip" exception in the action procedures:

              if "enter(len,o,c)" raises "Skip", "EnumPaths" will omit the 
              entire branch of the path tree rooted at current path, 
              and call "exit(len,o,c)" right away, as if "o" had no outgoing arcs;

              if "push(len,o,c,i,r,d)" raises "Skip", "EnumPaths" will
              call "pop(len,o,c,i,r,d)" right away;

              if "pop(len,o,c,i,r,d)" raises "Skip", "EnumPaths" ignores any remaining
              arcs out of "o", and calls "exit(len,o,c)" right away.

          (It is a checked error for "exit" to raise "Skip".)

          If any action is not specified, "EnumPaths" will provide a
          a trivial default action that does nothing (and raises no exceptions). *)

      (***************************************************************************)
      (* INPUT/OUTPUT                                                            *)
      (***************************************************************************)

      Dump(wr: Wr.T);
        (*
          Writes a description of the "DAG.T" to "wr", in a
          format that can be read back with "Load" below.  The output
          is in ASCII but not meant to be human-readable. *)

    END;

TYPE
  NodeAction = PROCEDURE(len: NAT; s: Node; cass: Symbol)
    RAISES {Abort, Skip};

  ArcAction = PROCEDURE(
      len: NAT; 
      org: Node; class: Symbol; 
      i: NAT; label: Symbol; dest: Node;
    ) RAISES {Abort, Skip};

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
