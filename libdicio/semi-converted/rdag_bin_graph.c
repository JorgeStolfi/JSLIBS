/* Last edited on 1999-06-05 20:50:50 by stolfi */

#ifndef _H
#define _H


/* See the copyright and disclaimer note at the end of this file. */

IMPORT Rd, Wr, Thread;
IMPORT FileFmt, FGet, FPut, NGet, NPut;
IMPORT BinGraphF;
FROM Basics IMPORT NAT;
FROM Basics IMPORT Full, Skip, Abort;
FROM BinGraphF IMPORT Entry, Entries; 

REVEAL
  T == BinGraphF.Private BRANDED OBJECT
    OVERRIDES
      /* Inherited from "BinGraph.Public" */
      NNodes    = NNodes;
      Get       = Get;
      EnumPaths = EnumPaths;
      Dump      = Dump;
      Load      = Load;

      /* Inherited from "BinGraphF.Private" */
      NAlloc    = NAlloc;
      Alloc     = Alloc;
      Set       = Set;
      Add       = Add;
      Clear     = Clear;
      Copy      = Copy;
      Crunch    = Crunch;
    ;};

PROCEDURE NNodes(g: T): Node ==
  {
    return g.nNodes
  ;} NNodes;

PROCEDURE Get(g: T; s: Node): NodeData ==
  {
    assert(s < g.nNodes );
    with (es == g.e[s]){
      return NodeData{
        lSym = es.lSym, lNode = es.lNode, 
        rSym = es.rSym, rNode = es.rNode
      }
    ;}
  ;} Get;

PROCEDURE NAlloc(g: T): NAT ==
  {
    if ((g.e == NULL)){ return 0 }else{ return NUMBER(g.e^) ;};
  ;} NAlloc;

PROCEDURE Alloc(g: T; maxNodes: NAT): T ==
  {
    assert(maxNodes >= g.nNodes );
    /* See if we already have the right amount of storage: */
    with (nAlloc == g.NAlloc()){
      if ((nAlloc >= maxNodes)  AND  AND  (nAlloc - maxNodes < maxNodes)){ return g ;};
      g.e = NULL;
      if ((maxNodes == 0)){ return g ;};
    ;};
    
    /* Needs to (re)allocate: */
    with (
      e == g.e,
      eNew == NEW(REF Entries, maxNodes)
   ){
      if ((g.nNodes > 0)){
        SUBARRAY(eNew^, 0, g.nNodes) = e^
      ;};
      e = eNew
    ;};
    return g
  ;} Alloc;

PROCEDURE Set(g: T; s: Node; lSym, rSym: Symbol; lNode, rNode: Node) ==
  {
    assert(s < g.nNodes );
    assert(lNode < g.nNodes );
    assert(rNode < g.nNodes );
    with (e == g.e^){
      with (
        es == e[s],
        en == Entry{
          es.lSym = lSym,
          es.lNode = lNode,
          es.rSym = rSym,
          es.rNode = rNode
        }
     ){
        if ((en!=es)){ 
          es = en;
          INC(g.epoch)
        ;}
      ;}
    ;}
  ;} Set;
  
PROCEDURE Add(g: T; lSym, rSym: Symbol; lNode, rNode: Node): Node RAISES {Full} ==
  {
    assert(s <= g.nNodes );
    if ((g.e == NULL) || (s > LAST(g.e^))){ RAISE Full ;};
    with (
      e == g.e^,
      s == g.nNodes + 0
   ){
      INC(g.nNodes);
      assert(lNode < g.nNodes );
      assert(rNode < g.nNodes );
      with (es == e[s],
        en == Entry{
          es.lSym = lSym,
          es.lNode = lNode,
          es.rSym = rSym,
          es.rNode = rNode
        }
     ){
        es = en
      ;};
      return s
    ;}
  ;} Add;
  
PROCEDURE Clear(g: T): T ==
  {
    if ((g.nNodes > 0)){
      g.nNodes = 0;
      g.e = NULL;
      INC(g.epoch)
    ;};
    return g
  ;} Clear;

PROCEDURE Copy(g: T; from: T; s: Node; VAR map: ARRAY OF Node): Node RAISES {Full} ==

  PROCEDURE DoCopy(u: Node): Node RAISES {Full} ==
    /*
      Does Copy starting at "u", a generic sucessor of "s": copies
      "u.lNode", then copies "u.rNode", and finally creates the copy of "u". */

    {
      if ((map[u] == NoNode)){
        map[u] = g.nNodes;
        with (
          eu == from.e[u],
          lNode == DoCopy(eu.lNode),
          rNode == DoCopy(eu.rNode)
       ){
          g.Set(map[u], eu.lSym, eu.rSym, lNode, rNode);
        ;}
      ;};
      return map[u]
    ;} DoCopy;

  {
    return DoCopy(s)
  ;} Copy;

PROCEDURE Crunch(
    T g; 
    READONLY root: ARRAY OF Node; 
    VAR /*OUT*/ map: ARRAY OF Node;
  ): T ==
  NAT *nOld;
      NAT nNew;
  {
    /*
      Crunch is done in in four passes:

      First pass: unmark all nodes ("map[s] == NoNode").

      Second pass: recursively mark ("map[s] == s") all 
      nodes reachable from "root".

      Third pass: scan all entries from 1 up, and move each marked
      node "s" to the lowest unused position, updating its "l" and
      "r" pointers, and saving that position in "map[s]".  Also
      update the root pointers.

      Fourth pass: rebuild the hash table for the compacted nodes.
    */

    if ((g.e == NULL)){ g.e = NEW(REF Entries, 0) ;};
    with (e == g.e^){

      /* Unmark all nodes */
      for (i = 0 TO g.nNodes-1){ map[i] = NoNode ;};

      /* Mark reachable nodes, and set "nOld > t" for all copied nodes "t" : */
      /* (convention: a node "t" is marked iff "map[t]==t".) */
      
      nOld = 0;
      
      PROCEDURE DoMark(t: Node) ==
        {
          while (map[t] == NoNode){
            if ((t >= nOld)){ nOld = t+1 ;};
            map[t] = t;
            with (lNode == e[t].lNode, rNode == e[t].rNode){
              if ((lNode < rNode)){
                DoMark(lNode); t = rNode
              }else{
                DoMark(rNode); t = lNode
              ;}
            ;}
          ;}
        ;} DoMark;
        
      {
        for (k = 0 TO LAST(root)){ DoMark(root[k]) ;}
      ;};
      
      /* Compress nodes virtually, updating "map": */
      nNew = 0;
      for (u = 0 TO nOld-1){
        if ((map[u]!=NoNode)){
          map[u] = nNew; INC(nNew)
        ;}
      ;};
      
      /* Now copy the node data, updating links: */
      for (u = 0 TO nOld-1){
        if ((map[u]!=NoNode)){
          with (
            eu == e[u], 
            lNode == map[eu.lNode], 
            rNode == map[eu.rNode],
            v == map[u], 
            ev == e[v]
         ){
            assert(lNode!=NoNode );
            assert(rNode!=NoNode );
            ev.lSym = eu.lSym;
            ev.rSym = eu.rSym;
            ev.lNode = lNode;
            ev.rNode = rNode
          ;};
        ;}
      ;};
      if ((nNew < g.nNodes)){
        INC(g.epoch);
        g.nNodes = nNew
      ;}
    ;};
    return g
  ;} Crunch;

PROCEDURE EnumPaths(
    T g;
    Node s;
    enter: NodeAction = NULL;
    push, pop: LinkAction = NULL;
    exit: NodeAction = NULL;
  ) RAISES {Abort} ==

  VAR len: NAT = 0;

  PROCEDURE DoEnumPaths(o: Node) RAISES {Abort} ==
    /*
      Does "EnumPaths" starting at "o", a generic sucessor of "s". */
    Node *t; x: Symbol;
    {
      with (eo == g.e[o]){
        TRY
          if ((enter!=NULL)){ enter(len, o, eo.lSym, eo.rSym) ;};
          for (b = FALSE TO TRUE){
            if ((NOT b)){ 
              t = eo.lNode; x = eo.lSym
            }else{
              t = eo.rNode; x = eo.lSym
            ;};
            TRY
              if ((push!=NULL)){ push(len, o, b, x, t) ;};
              INC(len);
              DoEnumPaths(t);
              DEC(len)
            EXCEPT 
              Skip ==> /*OK*/
            ;};
            if ((pop!=NULL)){ pop(len, o, b, x, t) ;};
          ;};
        EXCEPT
          Skip ==> /*OK*/
        ;};
        if ((exit!=NULL)){ 
          <*FATAL Skip ); { exit(len, o, eo.lSym, eo.rSym) ;}
        ;}
      ;}
    ;} DoEnumPaths;

  {
    assert(s < g.nNodes );
    assert(enter!=NULL) || (push!=NULL );
    DoEnumPaths(s)
  ;} EnumPaths;

CONST
  FileType == "BinGraph.T";
  FileVersion == "97-01-21";

  OldFileType     == "DAG.Dump";
  OldFileVersion  == "91-12-16";

PROCEDURE Dump(g: T; wr: Wr.T) ==
  <* FATAL Wr.Failure, Thread.Alerted );
  {
    FileFmt.WriteHeader(wr, FileType, FileVersion);
    FileFmt.WriteComment(wr, g.comment, '|');
    DumpBody(g, wr);
    FileFmt.WriteFooter(wr, FileType);
    fflush(wr);
  ;} Dump;

PROCEDURE DumpBody(g: T; wr: Wr.T) ==
  <* FATAL Wr.Failure, Thread.Alerted );
  {
    NPut.Int(wr, "nNodes", g.nNodes);
    with (e == g.e^){
      for (i = 0 TO g.nNodes-1){
        with (ei == e[i]){
          FPut.Int(wr, i);        FPut.Space(wr, 1);
          FPut.Int(wr, ei.lSym);  FPut.Space(wr, 1);
          FPut.Int(wr, ei.lNode); FPut.Space(wr, 1);
          FPut.Int(wr, ei.rSym);  FPut.Space(wr, 1);
          FPut.Int(wr, ei.rNode); FPut.EOL(wr);
        ;};
      ;};
    ;};
    fflush(wr);
  ;} DumpBody;

PROCEDURE Load(g: T; rd: Rd.T): T ==
  <* FATAL Rd.Failure, Rd.EndOfFile, Thread.Alerted );
  {
    with (hdr == Rd.GetLine(rd)){
      g.comment = FileFmt.ReadComment(rd, '|');
      if ((Text.Equal(hdr, FileFmt.MakeHeader(FileType, FileVersion)))){
        /* Official format. */
        LoadBody(g, rd);
        FileFmt.ReadFooter(rd, FileType)
      }else if ((Text.Equal(hdr, OldFileFmt.MakeHeader(OldFileType, OldFileVersion)))){
        /* Old "Dag.Dump" format. */
        LoadOldDAGBody(g, rd);
        OldFileFmt.ReadFooter(rd, OldFileType);
      ;};
    ;};
    return g
  ;} Load;

PROCEDURE LoadBody(g: T; rd: Rd.T) ==
  {
    with (nNodes == NGet.Int(rd, "nNodes")){
      /* Discard old contents and reallocate if necessary: */
      if ((g.nNodes > 0)){
        g.nNodes = 0;
        INC(g.epoch)
      ;};
      EVAL g.Alloc(nNodes);
      /* Read nodes: */
      with (e == g.e^){
        for (i = 0 TO nNodes-1){
          with (
            s == FGet.Int(rd),
            lSym  == FGet.Int(rd),
            lNode == FGet.Int(rd),
            rSym  == FGet.Int(rd),
            rNode == FGet.Int(rd)
         ){
            assert(s == i );
            assert(lNode < nNodes );
            assert(rNode < nNodes );
            e[s] = Entry{
              lSym  = lSym,
              lNode = lNode,
              rSym  = rSym,
              rNode = rNode
            };
          ;};
          FGet.EOL(rd);
        ;};
      ;};
      g.nNodes = nNodes;
    ;}
  ;} LoadBody;

PROCEDURE LoadOldDAGBody(g: T; rd: Rd.T) ==
  {
    with (nNodes == NGet.Int(rd, "max state") + 1){
      /* Discard old contents and reallocate if necessary: */
      if ((g.nNodes > 0)){
        g.nNodes = 0;
        INC(g.epoch)
      ;};
      EVAL g.Alloc(nNodes);
      with (e == g.e^){
        /* Create the "NullState": */
        e[0] = Entry{lSym = 0, lNode = 0, rSym = 0, rNode = 0};
        /* Read proper states: */
        for (i = 1 TO nNodes-1){
          with (
            s == FGet.Int(rd),
            rSym == FGet.Int(rd),
            lSym == FGet.Int(rd),
            rNode == FGet.Int(rd),
            lNode == FGet.Int(rd)
         ){
            assert(s == i );
            assert(lNode < nNodes );
            assert(rNode < nNodes );
            e[s] = Entry{
              lSym  = lSym,
              lNode = lNode,
              rSym  = rSym,
              rNode = rNode
            }
          ;};
          FGet.EOL(rd)
        ;}
      ;};
      g.nNodes = nNodes;
    ;}
  ;} LoadOldDAGBody;

{
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
