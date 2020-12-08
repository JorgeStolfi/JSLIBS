INTERFACE Code;

IMPORT Rd, Word, Wr;

TYPE
  TypeSize = [0 .. Word.Size];
  T       <: Public;
  Public   = OBJECT
    METHODS
      Fin(VAR t: TEXT);
      InitRead(rd: Rd.T);
      InitWrite(wr: Wr.T);
      Read(size: TypeSize): CARDINAL RAISES {Rd.EndOfFile};
      Write(size: TypeSize; n: CARDINAL)
    END;
  PROCEDURE BitSize(n: CARDINAL): CARDINAL;
  PROCEDURE BooleanSize(): CARDINAL;
  PROCEDURE NewCode(): T;
  PROCEDURE ViewBits(rd: Rd.T; wr: Wr.T);
  PROCEDURE WordSize(): CARDINAL;
  PROCEDURE WordSizeSize(): CARDINAL;
END Code.
