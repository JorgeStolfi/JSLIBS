/****************************************************************************/
/* (C) Copyright 1994 Universidade Estadual de Campinas (UNICAMP)           */
/*                    Campinas, SP, Brazil                                  */
/*                                                                          */
/* Authors:                                                                 */
/*                                                                          */
/*   Tomasz Kowaltowski  - CS Dept, UNICAMP <tomasz@dcc.unicamp.br>         */
/*   Claudio L. Lucchesi - CS Dept, UNICAMP <lucchesi@dcc.unicamp.br>       */
/*   Jorge Stolfi        - CS Dept, UNICAMP <stolfi@dcc.unicamp.br>         */
/*                                                                          */
/* This file can be freely distributed, modified, and used for any          */
/*   non-commercial purpose, provided that this copyright and authorship    */
/*   notice be included in any copy or derived version of this file.        */
/*                                                                          */
/* DISCLAIMER: This software is offered ``as is'', without any guarantee    */
/*   as to fitness for any particular purpose.  Neither the copyright       */
/*   holder nor the authors or their employers can be held responsible for  */
/*   any damages that may result from its use.                              */
/****************************************************************************/
#ifndef _H
#define _H


IMPORT Rd, Code;

EXCEPTION OverFlow;

TYPE
  FrequencyTable == REF ARRAY OF unsigned;
  T              <: Public;
  Public          == OBJECT
    METHODS
      Compress      (c: CHAR; cd: Code.T);
      Uncompress    (cd: Code.T): CHAR RAISES {Rd.EndOfFile};
      Dump          (cd: Code.T);
    ;};
  PROCEDURE New         (TblSiz: unsigned = NUMBER(CHAR)): T;
  PROCEDURE Tbl2Huffman (f: FrequencyTable): T RAISES {OverFlow};
  PROCEDURE Load        (cd: Code.T; TblSiz: unsigned = NUMBER(CHAR)): T 
    RAISES {Rd.EndOfFile};
  
;} Huffman.
