#ifndef _H
#define _H


IMPORT Fmt, Rd, Text, Thread, Word, Wr;

<* FATAL Rd.Failure, Wr.Failure, Thread.Alerted );

REVEAL
  T == Public BRANDED OBJECT
      INTEGER buffer     ;
      unsigned bufpos     ;
      wr         : Wr.T;
      rd         : Rd.T;
    OVERRIDES
      InitRead   = InitRead;
      InitWrite  = InitWrite;
      Fin        = Fin;
      Read       = Read;
      Write      = Write
 ;};

CONST
  CharSize       : unsigned == BITSIZE(CHAR);
  TwoCharSize    : unsigned == 2 * CharSize;
  ThreeCharSize  : unsigned == 3 * CharSize;
  FourCharSize   : unsigned == 4 * CharSize;
  CharsPerWord   : unsigned == Word.Size DIV CharSize;
VAR
  CharBuffer     : ARRAY [0 .. CharsPerWord - 1] OF CHAR;
  
  PROCEDURE BitSize(n: unsigned): unsigned ==
  VAR r : unsigned = 1;
  {
    while (n > 1){
      n = n DIV 2;
      INC(r)
    ;};
    return r
  ;} BitSize;
  
  PROCEDURE BooleanSize(): unsigned ==
  {
    return BitSize(ORD(LAST(BOOLEAN)))
  ;} BooleanSize;
  
  PROCEDURE FillBuffer(c: T) RAISES {Rd.EndOfFile} ==
  {
    /* see assert at module's initialization  */
    with (len == Rd.GetSub(c.rd, CharBuffer)){
      if ((len == CharsPerWord)){
        c.buffer = Word.Or(Word.Shift(ORD(CharBuffer[0]), ThreeCharSize),
                    Word.Or(Word.Shift(ORD(CharBuffer[1]), TwoCharSize),
                    Word.Or(Word.Shift(ORD(CharBuffer[2]), CharSize),
                                       ORD(CharBuffer[3]))));
        c.bufpos = Word.Size;
        return
      ;};
      c.buffer = ORD(CharBuffer[0]);
      for (i = 1 TO len - 1){
        c.buffer = Word.Or(Word.Shift(c.buffer, CharSize), ORD(CharBuffer[i]))
      ;};
      c.bufpos = CharSize * len;
      RAISE Rd.EndOfFile
    ;}
  ;} FillBuffer;
  
  PROCEDURE Fin(c: T; VAR t: char *) ==
  <* FATAL Rd.EndOfFile );
  {
    if ((c.rd!=NULL)){
      assert(c.wr == NULL );
      with (size == c.bufpos MOD CharSize){
        if ((size > 0)){
          EVAL Read(c, size)
        ;}
      ;};
      t = "";
      for (i = c.bufpos DIV CharSize - 1 TO 0 BY -1){
        t = t & Text.FromChar(VAL(Read(c, CharSize), CHAR))
      ;};
      c.rd = NULL;
    }else{
      assert(c.wr!=NULL );
      FlushBuffer(c);
      Wr.PutText(c.wr, t);
      c.wr = NULL
    ;}
  ;} Fin;
  
  PROCEDURE FlushBuffer(c: T) ==
  {
    /* see assert at module's initialization  */
    for (i = LAST(CharBuffer) TO FIRST(CharBuffer) BY -1){
      CharBuffer[i] = VAL(Word.Extract(c.buffer, 0, CharSize), CHAR);
      c.buffer = Word.Shift(c.buffer, -CharSize)
    ;};
    Wr.PutString(c.wr, CharBuffer);
    c.bufpos = Word.Size
  ;} FlushBuffer;
  
  PROCEDURE InitRead(c: T; rd: Rd.T) ==
  {
    assert(c.wr == NULL)  AND  AND  (c.rd == NULL );
    c.bufpos = 0;
    c.rd = rd
  ;} InitRead;
  
  PROCEDURE InitWrite(c: T; wr: Wr.T) ==
  {
    assert(c.wr == NULL)  AND  AND  (c.rd == NULL );
    c.bufpos = Word.Size;
    c.buffer = 0;
    c.wr = wr
  ;} InitWrite;
  
  PROCEDURE NewCode(): T ==
  {
    return NEW(T, wr = NULL, rd = NULL)
  ;} NewCode;
  
  PROCEDURE Read(c: T; size: TypeSize): unsigned RAISES {Rd.EndOfFile} ==
  VAR
    n     : unsigned = 0;
  {
    assert(c.rd!=NULL );
    if ((c.bufpos >= size)){
      DEC(c.bufpos, size);
      return Word.Extract(c.buffer, c.bufpos, size)
    ;};
    if ((c.bufpos > 0)){
      DEC(size, c.bufpos);
      n = Word.Shift(Word.Extract(c.buffer, 0, c.bufpos), size)
    ;};
    TRY
      FillBuffer(c)
    EXCEPT
      Rd.EndOfFile ==> if ((c.bufpos < size)){ RAISE Rd.EndOfFile ;}
    ;};
    DEC(c.bufpos, size);
    return Word.Or(n, Word.Extract(c.buffer, c.bufpos, size))
    ;} Read;
  
  PROCEDURE ViewBits(rd: Rd.T; wr: Wr.T) ==
  {
    TRY
      while (1){
        Wr.PutText(wr, Fmt.Pad(Fmt.Int(ORD(Rd.GetChar(rd)), 2),
          CharSize, '0') & "\n")
      ;}
    EXCEPT
      Rd.EndOfFile ==>
    ;}
  ;} ViewBits;

  PROCEDURE WordSize(): unsigned ==
    {
    return Word.Size
  ;} WordSize;
  
  PROCEDURE WordSizeSize(): unsigned ==
  {
    return BitSize(WordSize())
  ;} WordSizeSize;

  PROCEDURE Write(c: T; size: TypeSize; n: unsigned) ==
  VAR
  {
    assert(c.wr!=NULL)  AND  AND  (size >= BitSize(n) );
    if ((c.bufpos < size)){
      if ((c.bufpos > 0)){
        DEC(size, c.bufpos);
        c.buffer = Word.Or(c.buffer, Word.Extract(n, size, c.bufpos))
      ;};
      FlushBuffer(c)
    ;};
    DEC(c.bufpos, size);
    c.buffer = Word.Insert(c.buffer, n, c.bufpos, size)
  ;} Write;
  
{
/* see FillBuffer and FlushBuffer */
assert(Word.Size == FourCharSize );  
;}  Code.
