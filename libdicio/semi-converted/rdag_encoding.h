#ifndef _H
#define _H


/* Abstract mappings between Strings and ARRAY OF CHAR. */
/* See the copyright and disclaimer note at the end of this file. */

IMPORT Wr, Rd, Thread;
FROM Basics IMPORT Symbol, String, Done;
FROM Basics IMPORT NAT, CHARS;

EXCEPTION
  BadChar(CHAR);     /* Invalid character. */
  BadLetter(Symbol); /* Invalid Symbol. */
  BadChars(char *);    /* CHARS is untranslatable to String (arg tells why). */
  BadString(char *);   /* String is untranslatable to CHARS (arg tells why). */

TYPE
  T == OBJECT METHODS

      /* OUTPUT: STRING TO CHARS */

      StringToChars(
        String READONLY w; 
        VAR /*IO*/ c: REF CHARS; 
        VAR /*OUT*/ len: NAT;
      ) RAISES {BadString};
      /*
        Converts a given string "w" to an array of CHAR.
        The result is returned in SUBARRAY(c^, 0, len).
        Will automatically reallocate "c" if it is NULL or too short. */

      StringToText(READONLY w: String): char *RAISES {BadString};
      /*
        Converts a String "w" to a char *.  Equivalent to 
        |   StringToChars(w, c, len); 
        |   return Text.FromChars(SUBARRAY(c^, 0, len))
        but perhaps more efficient. */

      PrintString(
        wr: Wr.T; 
        String READONLY w;
      ) RAISES {BadString, Wr.Failure, Thread.Alerted};
      /*
        Equivalent to Wr.PutText(wr, StringToText(w)),
        but perhaps more efficient. */

      PrintLetter(
        wr: Wr.T;
        symbol: Symbol;
      ) RAISES {Wr.Failure, Thread.Alerted};
      /*
        Prints the given single "symbol" to "wr".

        This method is generally used when printing strings of symbols
        that may not be complete words: for instance, when printing
        isolated transitions of an automaton, or the suffixes and
        prefixes of a state, etc.

        Therefore, the PrintLetter method of any subtype should be
        able to handle any Symbol code, and print it in a legible and
        preferably unambiguous format.  

        It is not necessary that PrintString(wr, w) be equivalent to
        applying PrintLetter(wr, x) for each symbol of the string "w";
        indeed, with many reasonable encodings, this may be just
        impossible. */

      /* INPUT: CHARS TO STRING */

      /*
        All the procedures below return the result IN SUBARRAY(w^, 0, len).
        They will automatically reallocate "w" if the given one is
        NULL or too short. */

      CharsToString(
        CHARS READONLY c;
        VAR /*IO*/ w: REF String;
        VAR /*OUT*/ len: NAT;
      ) RAISES {BadChars};
      /*
        Converts a given array of CHAR "c" to a String. */

      TextToString(
        char *t;
        VAR /*IO*/ w: REF String;
        VAR /*OUT*/ len: NAT;
      ) RAISES {BadChars};
      /*
        Converts a char *"t" to a String */

      ReadString(
        rd: Rd.T;
        VAR /*IO*/ w: REF String;
        VAR /*OUT*/ len: NAT
      ) RAISES {BadChars, Done, Rd.Failure, Thread.Alerted};
      /*
        Reads a String from the given "rd". If there are no more strings
        in "rd", raises "Done".  

        If the next string is invalid, raises BadChars(msg)
        and leaves the reader positioned right AFTER the character
        that triggered the error. 

        The definition of a valid string and the criteria for deciding
        where a string ends and when there are no more strings
        are instance-specific. */

    ;};

;} Encoding.

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
