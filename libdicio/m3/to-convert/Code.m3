MODULE  Code;

IMPORT Fmt, Rd, Text, Thread, Word, Wr;

<* FATAL Rd.Failure, Wr.Failure, Thread.Alerted *>

REVEAL
  T = Public BRANDED OBJECT
      buffer     : INTEGER;
      bufpos     : CARDINAL;
      wr         : Wr.T;
      rd         : Rd.T;
    OVERRIDES
      InitRead   := InitRead;
      InitWrite  := InitWrite;
      Fin        := Fin;
      Read       := Read;
      Write      := Write
 END;

CONST
  CharSize       : CARDINAL = BITSIZE(CHAR);
  TwoCharSize    : CARDINAL = 2 * CharSize;
  ThreeCharSize  : CARDINAL = 3 * CharSize;
  FourCharSize   : CARDINAL = 4 * CharSize;
  CharsPerWord   : CARDINAL = Word.Size DIV CharSize;
VAR
  CharBuffer     : ARRAY [0 .. CharsPerWord - 1] OF CHAR;
  
  PROCEDURE BitSize(n: CARDINAL): CARDINAL =
  VAR r : CARDINAL := 1;
  BEGIN
    WHILE n > 1 DO
      n := n DIV 2;
      INC(r)
    END;
    RETURN r
  END BitSize;
  
  PROCEDURE BooleanSize(): CARDINAL =
  BEGIN
    RETURN BitSize(ORD(LAST(BOOLEAN)))
  END BooleanSize;
  
  PROCEDURE FillBuffer(c: T) RAISES {Rd.EndOfFile} =
  BEGIN
    (* see assert at module's initialization  *)
    WITH len = Rd.GetSub(c.rd, CharBuffer) DO
      IF len = CharsPerWord THEN
        c.buffer := Word.Or(Word.Shift(ORD(CharBuffer[0]), ThreeCharSize),
                    Word.Or(Word.Shift(ORD(CharBuffer[1]), TwoCharSize),
                    Word.Or(Word.Shift(ORD(CharBuffer[2]), CharSize),
                                       ORD(CharBuffer[3]))));
        c.bufpos := Word.Size;
        RETURN
      END;
      c.buffer := ORD(CharBuffer[0]);
      FOR i := 1 TO len - 1 DO
        c.buffer := Word.Or(Word.Shift(c.buffer, CharSize), ORD(CharBuffer[i]))
      END;
      c.bufpos := CharSize * len;
      RAISE Rd.EndOfFile
    END
  END FillBuffer;
  
  PROCEDURE Fin(c: T; VAR t: TEXT) =
  <* FATAL Rd.EndOfFile *>
  BEGIN
    IF c.rd # NIL THEN
      <* ASSERT c.wr = NIL *>
      WITH size = c.bufpos MOD CharSize DO
        IF size > 0 THEN
          EVAL Read(c, size)
        END
      END;
      t := "";
      FOR i := c.bufpos DIV CharSize - 1 TO 0 BY -1 DO
        t := t & Text.FromChar(VAL(Read(c, CharSize), CHAR))
      END;
      c.rd := NIL;
    ELSE
      <* ASSERT c.wr # NIL *>
      FlushBuffer(c);
      Wr.PutText(c.wr, t);
      c.wr := NIL
    END
  END Fin;
  
  PROCEDURE FlushBuffer(c: T) =
  BEGIN
    (* see assert at module's initialization  *)
    FOR i := LAST(CharBuffer) TO FIRST(CharBuffer) BY -1 DO
      CharBuffer[i] := VAL(Word.Extract(c.buffer, 0, CharSize), CHAR);
      c.buffer := Word.Shift(c.buffer, -CharSize)
    END;
    Wr.PutString(c.wr, CharBuffer);
    c.bufpos := Word.Size
  END FlushBuffer;
  
  PROCEDURE InitRead(c: T; rd: Rd.T) =
  BEGIN
    <* ASSERT c.wr = NIL AND c.rd = NIL *>
    c.bufpos := 0;
    c.rd := rd
  END InitRead;
  
  PROCEDURE InitWrite(c: T; wr: Wr.T) =
  BEGIN
    <* ASSERT c.wr = NIL AND c.rd = NIL *>
    c.bufpos := Word.Size;
    c.buffer := 0;
    c.wr := wr
  END InitWrite;
  
  PROCEDURE NewCode(): T =
  BEGIN
    RETURN NEW(T, wr := NIL, rd := NIL)
  END NewCode;
  
  PROCEDURE Read(c: T; size: TypeSize): CARDINAL RAISES {Rd.EndOfFile} =
  VAR
    n     : CARDINAL := 0;
  BEGIN
    <* ASSERT c.rd # NIL *>
    IF c.bufpos >= size THEN
      DEC(c.bufpos, size);
      RETURN Word.Extract(c.buffer, c.bufpos, size)
    END;
    IF c.bufpos > 0 THEN
      DEC(size, c.bufpos);
      n := Word.Shift(Word.Extract(c.buffer, 0, c.bufpos), size)
    END;
    TRY
      FillBuffer(c)
    EXCEPT
      Rd.EndOfFile => IF c.bufpos < size THEN RAISE Rd.EndOfFile END
    END;
    DEC(c.bufpos, size);
    RETURN Word.Or(n, Word.Extract(c.buffer, c.bufpos, size))
    END Read;
  
  PROCEDURE ViewBits(rd: Rd.T; wr: Wr.T) =
  BEGIN
    TRY
      LOOP
        Wr.PutText(wr, Fmt.Pad(Fmt.Int(ORD(Rd.GetChar(rd)), 2),
          CharSize, '0') & "\n")
      END
    EXCEPT
      Rd.EndOfFile =>
    END
  END ViewBits;

  PROCEDURE WordSize(): CARDINAL =
    BEGIN
    RETURN Word.Size
  END WordSize;
  
  PROCEDURE WordSizeSize(): CARDINAL =
  BEGIN
    RETURN BitSize(WordSize())
  END WordSizeSize;

  PROCEDURE Write(c: T; size: TypeSize; n: CARDINAL) =
  VAR
  BEGIN
    <* ASSERT c.wr # NIL AND size >= BitSize(n) *>
    IF c.bufpos < size THEN
      IF c.bufpos > 0 THEN
        DEC(size, c.bufpos);
        c.buffer := Word.Or(c.buffer, Word.Extract(n, size, c.bufpos))
      END;
      FlushBuffer(c)
    END;
    DEC(c.bufpos, size);
    c.buffer := Word.Insert(c.buffer, n, c.bufpos, size)
  END Write;
  
BEGIN
(* see FillBuffer and FlushBuffer *)
<* ASSERT Word.Size = FourCharSize *>  
END  Code.
