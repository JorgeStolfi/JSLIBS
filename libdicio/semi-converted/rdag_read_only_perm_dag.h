#ifndef _H
#define _H


IMPORT Basics, Rd;

TYPE
  State                 == unsigned;
  Letter                == Basics.Letter;
  ExtendedLetter        == [FIRST(Letter) - 1 .. LAST(Letter)];
CONST
  FinalBitChar          : CHAR     == '@';
  InexistentLetter      : INTEGER  == FIRST(ExtendedLetter);
  ComprHeader           : char *== "©KLS";
  ComprTrailer          : char *== "SLK©";
  DefaultCommentText    : char *== "no comment";

TYPE
  T <: Public;
  Public == OBJECT
    METHODS
      Accepts       (    s     : State; word: char *): BOOLEAN;
      ExpandState   (    s     : State;
                     BOOLEAN *final ;
                     VAR chars : REF ARRAY OF CHAR);
      IsFinal       (    s     : State): BOOLEAN;
      MakeTransition(    s     : State; c: CHAR): State;
      Root          (    level : unsigned = 0): State
    ;};
  PROCEDURE LoadCompr(rd: Rd.T): T;
;} ReadOnlyPermDAG.
