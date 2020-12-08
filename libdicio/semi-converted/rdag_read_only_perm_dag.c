#ifndef _H
#define _H


IMPORT Code, Huffman;
IMPORT Rd, FGet, Thread, Text;

CONST
  CharSize     == BITSIZE(CHAR);
VAR
  BooleanSize  : unsigned = Code.BooleanSize();
  WordSize     : unsigned = Code.WordSize();
  WordSizeSize : unsigned = Code.WordSizeSize() - 1; /* excludes 0 */
TYPE
  EntryRecord  == RECORD
      Letter rd   ;
      State dest ;
      rest : State
    ;};
REVEAL
    T == Public BRANDED OBJECT
      char *comment            ;
      unsigned fields             ;
      char *textperm           ;
      charperm           : REF ARRAY OF CHAR;
      chartolettertable  : ARRAY CHAR OF ExtendedLetter;
      root               : REF ARRAY OF State;
      dag                : REF ARRAY OF EntryRecord;
    OVERRIDES
      Accepts            = Accepts;
      ExpandState        = ExpandState;
      IsFinal            = IsFinal;
      MakeTransition     = MakeTransition;
      Root               = Root
    ;};

  PROCEDURE Accepts(permdag: T; s: State; word: char *): BOOLEAN ==
  {
    with (len == Text.Length(word),
         w   == word & Text.FromChar(FinalBitChar),
         dag == permdag.dag,
         t   == permdag.chartolettertable
        ){
      for (i = 0 TO len){
        with (c == t[Text.GetChar(w, i)]){
          while (1){
            if ((s == 0)){ return FALSE ;};
            with (e == dag[s]){
              if ((e.rd == c)){ s = e.dest; EXIT
              }else{ s = e.rest
              ;}
            ;}
          ;}
        ;}
      ;};
      return s == 0
    ;}
  ;} Accepts;
  
  PROCEDURE ExpandState(    permdag : T; 
                            State s       ;
                        BOOLEAN *final   ;
                        VAR chars   : REF ARRAY OF CHAR) ==
  VAR
    fullchars : ARRAY CHAR OF CHAR;
    nchars    : CHAR = FIRST(CHAR);
  {
    final = FALSE;
    with (dag  == permdag.dag,
         perm == permdag.charperm
        ){
      while (s > 0){
        with (e == dag[s],
             c   == perm[e.rd]
            ){
          if ((c == FinalBitChar)  AND  AND  (e.dest == 0)){
            final = TRUE
          }else{
            fullchars[nchars] = c;
            INC(nchars)
          ;};
          s = e.rest
        ;}
      ;}
    ;};
    with (n == ORD(nchars)){
      chars  = NEW(REF ARRAY OF CHAR, n);
      chars^ = SUBARRAY(fullchars, 0, n)
    ;}
  ;} ExpandState;
  
  PROCEDURE IsFinal(permdag: T; s: State): BOOLEAN ==
  VAR
    BOOLEAN final ;
    chars : REF ARRAY OF CHAR;
  {
    ExpandState(permdag, s, final, chars);
    return final
  ;} IsFinal;
  
  PROCEDURE LoadCompr(rd: Rd.T): T ==
  <* FATAL Rd.EndOfFile, Rd.Failure, Thread.Alerted );
  VAR
    permdag      : T = NEW(T);
    firststate   : State = 1;
    State laststate    ;
    State dest         ;
    State rest         ;
    
    code         : Code.T = Code.NewCode();
    Letter rdsize       ;
    unsigned msize        ;
    char *tail         ;
    h            : Huffman.T;
    unsigned rdVar        ;
  {
    FGet.Match(rd, ComprHeader); FGet.EOL(rd);
    code.InitRead(rd);
    with (fields == permdag.fields,
         com    == permdag.comment){
      fields = code.Read(WordSizeSize);
      if ((fields >= 1)){
        com = "";
        for (i = 1 TO code.Read(WordSize)){
          com = com & Text.FromChar(VAL(code.Read(CharSize), CHAR))
        ;}
      }else{
        com = DefaultCommentText
      ;}
    ;};
    rdsize = code.Read(CharSize);
    with (charperm   == permdag.charperm){
      charperm = NEW(REF ARRAY OF CHAR, 1 + code.Read(rdsize));
      for (c = FIRST(CHAR) TO LAST(CHAR)){
        with (codet == VAL(code.Read(BooleanSize), BOOLEAN),
             t     == permdag.chartolettertable[c]
            ){
          if ((codet)){
            t = code.Read(rdsize);
            charperm[t] = c
          }else{
            t = InexistentLetter
          ;}
        ;}
      ;};
      permdag.textperm = Text.FromChars(charperm^)
    ;};   
    msize = 1 + code.Read(WordSizeSize);
    with (root == permdag.root){
      root = NEW(REF ARRAY OF State, code.Read(WordSizeSize));
      for (i = 0 TO LAST(root^)){
        root[i] = code.Read(msize)
      ;}
    ;};

    if ((permdag.fields >= 4)){ 
      h = Huffman.Load(code, NUMBER(permdag.charperm^));
    ;};

    with (m     == code.Read(msize),
         nulrd == permdag.chartolettertable[FinalBitChar],
         dag   == permdag.dag 
        ){
      dag = NEW(REF ARRAY OF EntryRecord, 1 + m);
      laststate = MIN(m, 2);
      for (statesize = 0 TO msize){
        for (s = firststate TO laststate){
          if ((permdag.fields >= 4 
         )){ rdVar = ORD(h.Uncompress(code))
          }else{ rdVar = code.Read(rdsize)
          ;};
          if ((VAL(code.Read(BooleanSize), BOOLEAN))){
            dest = 1 + code.Read(statesize)
          }else if ((rdVar == nulrd)){
            dest = 0
          }else{
            dest = s - 1
          ;};
          if ((VAL(code.Read(BooleanSize), BOOLEAN))){
            rest = s - 1
          }else if ((VAL(code.Read(BooleanSize), BOOLEAN))){
            rest = 1 + code.Read(statesize)
          }else{
            rest = 0
          ;};
          dag[s] = EntryRecord{rd = rdVar, dest = dest, rest= rest}
        ;};
        if ((laststate == m)){
          code.Fin(tail);
          EXIT
        }else{
          firststate = 1 + laststate;
          laststate = MIN(m, 2 * laststate - 1)
        ;}
      ;}
    ;};
    if ((permdag.fields >= 2)){
      if ((permdag.fields >= 3)){
        tail = ""
      ;};
    ;};

    /* Check trailer: */
    TRY
      with (line == Rd.GetLine(rd)){ tail = tail & line ;}
    EXCEPT 
      Rd.EndOfFile ==> /*OK*/
    ;};
    assert(Text.Equal(tail, ComprTrailer) );

    return permdag
  ;} LoadCompr;
  
  PROCEDURE MakeTransition(permdag: T; s: State; c: CHAR): State ==
  {
    with (dag    == permdag.dag,
         letter == permdag.chartolettertable[c]
        ){
      while (s > 0){
        with (e == dag[s]){
          if ((e.rd == letter)){ return e.dest ;};
          s = e.rest
        ;}
      ;}
    ;};
    return 0
  ;} MakeTransition;
  
  PROCEDURE Root(permdag: T; level : unsigned = 0): State ==
  {
    with (root == permdag.root){
      assert(level <= LAST(root^) );
      return root[level]
    ;}
  ;} Root;
  
{
;} ReadOnlyPermDAG.
