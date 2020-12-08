#ifndef _H
#define _H


IMPORT Rd, Word, Wr;

TYPE
  TypeSize == [0 .. Word.Size];
  T       <: Public;
  Public   == OBJECT
    METHODS
      Fin(VAR t: char *);
      InitRead(rd: Rd.T);
      InitWrite(wr: Wr.T);
      Read(size: TypeSize): unsigned RAISES {Rd.EndOfFile};
      Write(size: TypeSize; n: unsigned)
    ;};
  PROCEDURE BitSize(n: unsigned): unsigned;
  PROCEDURE BooleanSize(): unsigned;
  PROCEDURE NewCode(): T;
  PROCEDURE ViewBits(rd: Rd.T; wr: Wr.T);
  PROCEDURE WordSize(): unsigned;
  PROCEDURE WordSizeSize(): unsigned;
;} Code.
