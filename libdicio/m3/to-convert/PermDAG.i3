(****************************************************************************)
(* (C) Copyright 1992 Universidade Estadual de Campinas (UNICAMP)           *)
(*                    Campinas, SP, Brazil                                  *)
(*                                                                          *)
(* Authors:                                                                 *)
(*                                                                          *)
(*   Tomasz Kowaltowski  - CS Dept, UNICAMP <tomasz@dcc.unicamp.br>         *)
(*   Claudio L. Lucchesi - CS Dept, UNICAMP <lucchesi@dcc.unicamp.br>       *)
(*   Jorge Stolfi        - CS Dept, UNICAMP <stolfi@dcc.unicamp.br>         *)
(*                                                                          *)
(* This file can be freely distributed, modified, and used for any          *)
(*   non-commercial purpose, provided that this copyright and authorship    *)
(*   notice be included in any copy or derived version of this file.        *)
(*                                                                          *)
(* DISCLAIMER: This software is offered ``as is'', without any guarantee    *)
(*   as to fitness for any particular purpose.  Neither the copyright       *)
(*   holder nor the authors or their employers can be held responsible for  *)
(*   any damages that may result from its use.                              *)
(****************************************************************************)
INTERFACE PermDAG;

IMPORT Rd, Wr;

IMPORT Huffman;

FROM Basics IMPORT Full;

FROM ReadOnlyPermDAG IMPORT State;

TYPE
  T       <: Public;
  Public  = OBJECT
    METHODS
      Accepts       (    s           : State;
                         word        : TEXT)
                                     : BOOLEAN;
      
      AddSub        (    add         : BOOLEAN;
                         word        : TEXT;
                         level       : CARDINAL := 0)
                                       RAISES {Full};
      
      Comment       (    comment     : TEXT);
      
      Conta         (    wr          : Wr.T);
      
      Crunch        ();
      
      Dump          (    wr          : Wr.T);
      
      DumpCompr     (    wr          : Wr.T) RAISES {Huffman.OverFlow};
      
      DumpFixedBin  (    wr          : Wr.T);
      
      ExpandState   (    s           : State;
                     VAR final       : BOOLEAN;
                     VAR chars       : REF ARRAY OF CHAR);
      Fold          (    backupname  : TEXT;
                         backuplevel : CARDINAL;
                         restore     : BOOLEAN)
                                       RAISES {Full, Huffman.OverFlow};
                                       
      IncRoots      ();
      
      IsFinal       (    s           : State)
                                     : BOOLEAN;
                                     
      MakeTransition(    s           : State;
                         c           : CHAR)
                                     : State;
                                     
      Prm2Red       (    wr          : Wr.T)
                                       RAISES {Full};
                                       
      Root          (    level       : CARDINAL := 0)
                                     : State;
                                     
      SortPrm       ()                 RAISES {Full};
      
      Spell         (    wr          : Wr.T;
                         level       : CARDINAL := 0);
                         
      Statistics    (    wr          : Wr.T);
      
      UnFold        ()                 RAISES {Full}
    END;
  
  PROCEDURE Load     (rd: Rd.T): T RAISES {Full};
  PROCEDURE LoadCompr(rd: Rd.T): T RAISES {Full};
  PROCEDURE New      ()        : T;        
  PROCEDURE Red2Prm  (rd: Rd.T): T RAISES {Full};
END PermDAG.
