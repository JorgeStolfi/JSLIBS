#ifndef _H
#define _H


/* Trivial chars/string encoding: symbol == ORD(char). */
/* See the copyright and disclaimer note at the end of this file. */

/*
  This is the trivial character encoding "symbol == ORD(char)".
  It includes only the printing ISO-8859-1 chars 
  (ranges ['\040'..'\176'] and ['\241'..'\377']),
  excluding therefore all control characters, TAB, newline, etc.. 
  The encoded "Symbol" is always in [32..126] or [161..255].

  The "PrintLetter" method uses the Modula-3 escaped-octal 
  convention '\nnn' for non-printing characters.

  This particular encoding supports two additional methods,
  "SymbolToChar" and "CharToLetter".  They perform the 
  trivial encoding, but raise "BadChar" or "BadLetter"
  for non-printing characters.

  The "ReadString" method eliminates trailing and leading blanks, but
  preserves embedded blanks. The "CharsToString" and "TextToString"
  methods preserve all blanks. 
*/

IMPORT Encoding;

FROM Basics IMPORT Symbol;
FROM Encoding IMPORT BadChar, BadLetter;

CONST
  ValidChars == SET OF CHAR{'\040'..'\176', '\241'..'\255'};
  ValidLetters == SET OF Symbol{32..126, 161..255};

TYPE
  T <: Public;

  Public == Encoding.T OBJECT METHODS

      SymbolToChar(let: Symbol): CHAR RAISES {BadLetter};
      /*
        Converts a symbol to a single character */

      CharToLetter(ch: CHAR): Symbol RAISES {BadChar};
      /*
        Converts a single character to a symbol */

    ;};

PROCEDURE Init(e: T);
PROCEDURE New(): T;

;} PlainEncoding.

/****************************************************************************/
/* (C) Copyright 1992 Universidade Estadual de Campinas (UNICAMP)           */
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
