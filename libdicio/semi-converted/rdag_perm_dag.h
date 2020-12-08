#ifndef acau_perm_dag_H
#define acau_perm_dag_H


IMPORT Rd, Wr;

IMPORT Huffman;

FROM Basics IMPORT Full;

FROM ReadOnlyPermDAG IMPORT State;

TYPE
  T       <: Public;
  Public  == OBJECT
    METHODS
      Accepts       (    s           : State;
                         word        : char *)
                                    BOOLEAN  ;
      
      AddSub        (    add         : BOOLEAN;
                         char *word        ;
                         level       : unsigned = 0)
                                       RAISES {Full};
      
      Comment       (    comment     : char *);
      
      Conta         (    wr          : Wr.T);
      
      Crunch        ();
      
      Dump          (    wr          : Wr.T);
      
      DumpCompr     (    wr          : Wr.T) RAISES {Huffman.OverFlow};
      
      DumpFixedBin  (    wr          : Wr.T);
      
      ExpandState   (    s           : State;
                     BOOLEAN *final       ;
                     VAR chars       : REF ARRAY OF CHAR);
      Fold          (    backupname  : char *;
                         unsigned backuplevel ;
                         restore     : BOOLEAN)
                                       RAISES {Full, Huffman.OverFlow};
                                       
      IncRoots      ();
      
      IsFinal       (    s           : State)
                                    BOOLEAN  ;
                                     
      MakeTransition(    s           : State;
                         c           : CHAR)
                                    State  ;
                                     
      Prm2Red       (    wr          : Wr.T)
                                       RAISES {Full};
                                       
      Root          (    level       : unsigned = 0)
                                    State  ;
                                     
      SortPrm       ()                 RAISES {Full};
      
      Spell         (    wr          : Wr.T;
                         level       : unsigned = 0);
                         
      Statistics    (    wr          : Wr.T);
      
      UnFold        ()                 RAISES {Full}
    ;};
  
  PROCEDURE Load     (rd: Rd.T): T RAISES {Full};
  PROCEDURE LoadCompr(rd: Rd.T): T RAISES {Full};
  PROCEDURE New      ()        : T;        
  PROCEDURE Red2Prm  (rd: Rd.T): T RAISES {Full};
;} PermDAG.
