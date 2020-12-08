/* Last edited on 1999-06-05 20:40:12 by stolfi */

#ifndef _H
#define _H


/* The "BinGraph" data structure for labeled binary directed graphs. */
/* See the copyright and disclaimer note at the end of this file. */

/*
  A "BinGraph.T" is a labeled binary directed graph.

  Each node "s" of a "BinGraph.T" has two outgoing links, going
  to nodes "LNode(s)" and "RNode(s)", and two labels "LSym(s), RSym(s)".
  
  This interface gives only opaque, read-only access to a "BinGraph.T".
  See the ``friends'' interface "BingGraphF" for implementation details and 
  operations that modify the graph.
*/

IMPORT Basics, Commented, Wr;
FROM Basics IMPORT NAT, BOOL;
FROM Basics IMPORT Skip, Abort;

TYPE
  Symbol == Basics.Symbol;

  Node == NAT;  /* Index of a node of a "BinGraph.T" */
    /*
      Nodes are represented as "NAT"s, rather than opaque "REF"s,
      so that clients can use them as indices in tables, etc.

      Beware that a "BinGraph.Node" does *NOT* know which "BinGraph.T" it
      belongs to; the same "Node" value identifies unrelated nodes in
      different "BinGraph.T"s.  */

  NodeData == RECORD  /* Data for a particular node: */
      lSym: Symbol;  /* The "L" label */
      rSym: Symbol;  /* The "R" label */
      Node lNode;   /* Destination of the "L" link */
      Node rNode;   /* Destination of the "R" link */
    ;};

CONST 
  NoNode == LAST(Node); /* A "NULL" value for nodes. */

TYPE
  T <: Public;
  Public == Commented.T OBJECT
    METHODS

      NNodes(): Node;
        /*
          Nodes of this "BinGraph.T" are numbered from 0 to "NNodes()-1". */

      Get(s: Node): NodeData;
        /*
          The labels "LSym(s), RSym(s)", and the descendants "LNode(s)" 
          and "RNode(s)" of node "s".  Requires "s < NNodes()". 
          Cost: $O(1)$ time, 0 space. */

      /*******************************************************************/
      /* ENUMERATION                                                     */
      /*******************************************************************/

      EnumPaths(
          Node s;
          enter: NodeAction = NULL;
          push, pop: LinkAction = NULL;
          exit: NodeAction = NULL;
        ) RAISES {Abort};
        /* 
          Enumerates a set of paths in the "BinGraph.T", starting
          from node "s", in depth-first order. The set of paths is
          defined by the client-provided procedures "enter", "push",
          "pop", and "exit".

          The enumeration algorithm can be described in terms of a conceptual
          "current path" that initially contains just the node "s",
          and grows or shrinks one link at a time.  EnumPaths calls

              "enter(len,o,l,r)" whenever the current path has "len" ars and
                 has just reached node "o", whose symbol labels are "l,r".

              "push(len,o,b,x,t)" whenever the current path has "len"
                 links, ends at node "o", and is being extended to
                 node "t" through an link labeled "x"; where 
                 "t == LNode(o)" and "x == LSym(o)" if "b == FALSE", or 
                 "t == RNode(o)" and "x == RSym(o)" if "b == TRUE";

              "pop(len,o,b,x,t)" whenever the current path has "len+1" links,
                 ends at node "t", and EnumPaths is about to remove its
                 last link, which was labeled "x", and is either 
                 the "L" or "R" link (depending on "b") out of node "o";

              "exit(len,o,l,r)" whenever the current path has "len" links
                 and is about to pull back from node "o".

          Any of these client actions can stop the enumeration by raising "Abort".

          Unless the enumeration is aborted, every call to "enter" will
          be followed eventually by a call to "exit", and every call to
          "push" will be followed eventually by a matching call to
          "pop".  The typical action pattern for a generic node "o"

          |   enter(len,o,LSym(o),RSym(o));
          |     push(len,o,FALSE,LSym(o),LNode(o));
          |       ...
          |     pop(len,o,FALSE,LSym(o),LNode(o));
          |     push(len,o,TRUE,RSym(o),RNode(o));
          |       ...
          |     pop(len,o,TRUE,RSym(o),RNode(o));
          |   exit(len,o,LSym(o),RSym(o))

          In particular, "EnumPaths" will call
          "enter(0,s,RSym(s),LSym(s))" at the very beginning, and (if
          not aborted) "exit(0,s,RSym(s),LSym(s))" at the very end.

          Note that the same node "o" may be "enter"ed and "exit"ed many, many times.

          The client can prune branches of the path tree by 
          raising the "Skip" exception inside the action procedures.
          Namely,

              if "enter(len,o,r,l)" raises "Skip", "EnumPaths" will
              call "exit(len,o,r,l)" right away, ignoring the links 
              out of "o";

              if "push(len,o,b,x,t)" raises "Skip", "EnumPaths" will
              call "pop(len,o,b,x,t)" right away, ignoring the link in
              question;

              if "pop(len,o,b,x,t)" raises "Skip", "EnumPaths" ignores
              any remaining links out of "o" are ignored, and calls
              "exit(len,o,r,l)" right away.

          (It is a checked error for "exit" to raise "Skip".)

          If any action is not specified, "EnumPaths" will provide a
          a trivial default action that does nothing (and raises no exceptions). */

      /***************************************************************************/
      /* INPUT/OUTPUT                                                            */
      /***************************************************************************/

      Dump(wr: Wr.T);
        /*
          Writes a description of the graph to "wr".  The output is ASCII 
          text but not meant to be human-readable. */

    ;};

TYPE
  NodeAction == PROCEDURE(len: NAT; s: Node; l, r: Symbol)
    RAISES {Abort, Skip};

  LinkAction == PROCEDURE(len: NAT; o: Node; b: BOOL; x: Symbol; t: Node)
    RAISES {Abort, Skip};

;} BinGraph.

/****************************************************************************/
/* (C) Copyright 1992 Universidade Estadual de Campinas (UNICAMP)           */
/*                    Campinas, SP, Brazil                                  */
/*                                                                          */
/* Authors:                                                                 */
/*                                                                          */
/*   Tomasz Kowaltowski  - CS Dept, UNICAMP <tomasz@dcc.unicamp.br>         */
/*   Claudio L. Lucchesi - CS Dept, UNICAMP <lucchesi@dcc.unicamp.br>       */
/*   Jorge Stolfi        - CS Dept, UNICAMP <stolfi@dcc.unicamp.br>         */
/*                                                                          */
/* This file can be freely distributed, modified, and used for any          */
/*   non-commercial purpose, provided that this copyright and authorship    */
/*   notice be included in any copy or derived version of this file.        */
/*                                                                          */
/* DISCLAIMER: This software is offered ``as is'', without any guarantee    */
/*   as to fitness for any particular purpose.  Neither the copyright       */
/*   holder nor the authors or their employers can be held responsible for  */
/*   any damages that may result from its use.                              */
/****************************************************************************/
