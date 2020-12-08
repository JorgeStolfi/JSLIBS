INTERFACE ReadOnlyPermDAG;

IMPORT Basics, Rd;

TYPE
  State                 = CARDINAL;
  Letter                = Basics.Letter;
  ExtendedLetter        = [FIRST(Letter) - 1 .. LAST(Letter)];
CONST
  FinalBitChar          : CHAR     = '@';
  InexistentLetter      : INTEGER  = FIRST(ExtendedLetter);
  ComprHeader           : TEXT = "©KLS";
  ComprTrailer          : TEXT = "SLK©";
  DefaultCommentText    : TEXT = "no comment";

TYPE
  T <: Public;
  Public = OBJECT
    METHODS
      Accepts       (    s     : State; word: TEXT): BOOLEAN;
      ExpandState   (    s     : State;
                     VAR final : BOOLEAN;
                     VAR chars : REF ARRAY OF CHAR);
      IsFinal       (    s     : State): BOOLEAN;
      MakeTransition(    s     : State; c: CHAR): State;
      Root          (    level : CARDINAL := 0): State
    END;
  PROCEDURE LoadCompr(rd: Rd.T): T;
END ReadOnlyPermDAG.
