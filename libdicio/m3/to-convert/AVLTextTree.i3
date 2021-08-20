INTERFACE AVLTextTree;

  FROM Basics IMPORT Abort;

  TYPE
  
    T <: AVLPublic;
    AVLTextTreeState <: REFANY;

    AVLPublic = OBJECT
        METHODS
          Insert(READONLY w: TEXT);
          Search(READONLY w: TEXT): BOOLEAN;
          Delete(READONLY w: TEXT);
          Size(): CARDINAL;
          Enum(Action: TextActionProc; reverse: BOOLEAN := FALSE) RAISES{Abort};
          MakeTransition(s: AVLTextTreeState; ch: CHAR): AVLTextTreeState;
          IsFinal(s:AVLTextTreeState): BOOLEAN;
          RootState(): AVLTextTreeState;
        END;

      TextActionProc = PROCEDURE (READONLY w: TEXT) RAISES {Abort} ;
      
  PROCEDURE New(): T;

  PROCEDURE IsNullState(s: AVLTextTreeState): BOOLEAN;
  
END AVLTextTree.
