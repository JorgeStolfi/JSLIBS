#ifndef _H
#define _H


IMPORT Rd, Code, Word;

TYPE
  Bit == [0..1];
  Tree == REF TreeNode;
  TreeNode       == RECORD
      s: ARRAY Bit OF Tree = ARRAY Bit OF Tree{NULL, NULL};
      CHAR c;
    ;};
  CompressionTableEntry == RECORD
      l: Code.TypeSize;
      n: Word.T;
    ;};
  CodedTreeRecord       == RECORD
      a: REF ARRAY OF BOOLEAN;
      l: REF ARRAY OF CHAR;
    ;};
REVEAL
  T == Public BRANDED OBJECT
      unsigned TableSize        ;
      unsigned CharSize         ;
      Tree Root             ;
      CompressionTable : REF ARRAY OF CompressionTableEntry;
      unsigned TreeSize         ;
      CodedTree        : CodedTreeRecord
    OVERRIDES
      Compress         = Compress;
      Uncompress       = Uncompress;
      Dump             = Dump
    ;};
  
VAR
  BooleanSize = Code.BooleanSize();
  
  PROCEDURE New(TblSiz: unsigned): T ==
  {
    with (h == NEW(T)){
      h.TableSize        = TblSiz;
      h.CharSize         = Code.BitSize(h.TableSize - 1);
      h.TreeSize         = 1 + 2 * h.TableSize;
      h.Root             = NULL;
      h.CompressionTable = NEW(REF ARRAY OF CompressionTableEntry, h.TableSize);
      h.CodedTree        = CodedTreeRecord
        {a = NEW(REF ARRAY OF BOOLEAN, h.TreeSize),
         l = NEW(REF ARRAY OF CHAR,  h.TableSize)};
      return h
    ;}
  ;} New;

  PROCEDURE BuildTable(h: T) RAISES {OverFlow} ==

    PROCEDURE AuxBuild(t: Tree; ct: CompressionTableEntry) RAISES 
      {OverFlow} ==
      {
        with (e == t.s[0]){
          if ((e == NULL)){
            h.CompressionTable[ORD(t.c)] = ct
          }else if ((ct.l == Word.Size)){
            RAISE OverFlow
          }else{
            INC(ct.l); ct.n = Word.Shift(ct.n, 1); AuxBuild(e, ct);
            ct.n = Word.Or(ct.n, 1); AuxBuild(t.s[1], ct)
          ;}
        ;};
      ;} AuxBuild;

  {
    AuxBuild(h.Root, CompressionTableEntry{0, 0});
  ;} BuildTable;

  PROCEDURE EncodeTree(h: T) RAISES {} ==
  VAR
    na, nl : unsigned = 0;
    
    PROCEDURE AuxEncode(t: Tree) ==
    {
      with (e == t.s[0],
           x == e == NULL 
          ){
        h.CodedTree.a[na]  = x; INC(na);
        if ((x)){
          h.CodedTree.l[nl] = t.c; INC(nl)
        }else{
          AuxEncode(e); AuxEncode(t.s[1])
        ;}
      ;}
    ;} AuxEncode;
    
  {
    AuxEncode(h.Root)
  ;} EncodeTree;

  PROCEDURE Tbl2Huffman (f: FrequencyTable): T RAISES {OverFlow} ==
  VAR
    h: T = New(NUMBER(f^));
    FreqTab: REF ARRAY OF unsigned;
    NodeTab: REF ARRAY OF Tree;

    PROCEDURE Insert(n: unsigned) RAISES {} ==
      INTEGER *j;
      {
        assert(n > 0 );
        VAR 
          fr = FreqTab[n];
          no = NodeTab[n];
        {
          j = 0;
          for (i = n TO 1 BY -1){
            with (fi == FreqTab[i-1]){
              if ((fi <= fr)){ 
                FreqTab[i] = fi;
                NodeTab[i] = NodeTab[i-1]
              }else{ 
                j = i; 
                EXIT
              ;}
            ;}
          ;};
          FreqTab[j] = fr;
          NodeTab[j] = no;
        ;};
      ;} Insert;

  {
    FreqTab = NEW(REF ARRAY OF unsigned, h.TableSize);
    NodeTab = NEW(REF ARRAY OF Tree,     h.TableSize);
    for (n = 0 TO LAST(f^)){
      FreqTab[n] = f[n];
      NodeTab[n] = NEW(Tree, c = VAL(n, CHAR));
    ;};
    for (n = 1 TO LAST(f^)){ Insert(n) ;};
    for (n = LAST(f^) TO 1 BY -1){
      Insert(n);
      INC(FreqTab[n-1], FreqTab[n]);
      with (p == NodeTab[n-1],
           x == NEW(Tree)){
        x^ = p^;
        p.s = ARRAY Bit OF Tree{x, NodeTab[n]}
      ;}
    ;};
    h.Root = NodeTab[0];
    BuildTable(h);
    EncodeTree(h);
    return h
  ;} Tbl2Huffman;
  
  PROCEDURE Load        (cd: Code.T; TblSiz: unsigned = NUMBER(CHAR)): T 
    RAISES {Rd.EndOfFile} ==
  VAR
    h: T = New(TblSiz);
    na, nl : unsigned = 0;

    PROCEDURE AuxLoad(): Tree ==
    VAR
      t: Tree = NEW(Tree);
      BOOLEAN tr ;
    {
      tr = h.CodedTree.a[na]; INC(na);
      if ((tr)){
        t.c = h.CodedTree.l[nl]; INC(nl)
      }else{
        t.s[0] = AuxLoad(); t.s[1] = AuxLoad()
      ;};
      return t
    ;} AuxLoad;
    
  {
    for (i = 0 TO h.TreeSize - 1){
      with (a == h.CodedTree.a){
        a[i] = VAL(cd.Read(BooleanSize), BOOLEAN)
      ;}
    ;};
    for (i = 0 TO h.TableSize - 1){
      with (l == h.CodedTree.l){
        l[i] = VAL(cd.Read(h.CharSize), CHAR);
      ;}
    ;};
    h.Root = AuxLoad();
    return h
  ;} Load;
  
  PROCEDURE Compress      (h: T; c: CHAR; cd: Code.T) ==
  {
    with (ct == h.CompressionTable[ORD(c)]){
      cd.Write(ct.l, ct.n)
    ;}
  ;} Compress;
  
  PROCEDURE Uncompress    (h: T; cd: Code.T): CHAR RAISES {Rd.EndOfFile} ==
  VAR
    t: Tree = h.Root;
  {
    while (t.s[0]!=NULL ){
      t = t.s[cd.Read(1)];
    ;};
    return t.c
  ;} Uncompress;
  
  PROCEDURE Dump          (h: T; cd: Code.T) ==
  {
    for (i = 0 TO h.TreeSize - 1){ cd.Write(BooleanSize, ORD(h.CodedTree.a[i]))
    ;}; 
    for (i = 0 TO h.TableSize - 1){ cd.Write(h.CharSize, ORD(h.CodedTree.l[i]))
    ;}
  ;} Dump;

{
;} Huffman.
