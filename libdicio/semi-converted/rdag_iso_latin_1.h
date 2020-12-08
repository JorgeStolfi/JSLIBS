#ifndef _H
#define _H


/* The ISO Latin-1 character set, and related tools. */
/* See the copyright and disclaimer note at the end of this file. */
/* Last edited on 1999-06-05 20:56:36 by stolfi */

/*
  This interface defines the semantics of CHAR accoridng to the
  ISO-8859 (ISO Latin-1) standard.  It also provides mapping tables
  that translate lower-case symbols into upper-case and the like.

  Based on "ISOChar.i3" by Jim Meehan and Henri Gouraud 
  (Copyright 1994, Digital Equipment Corp.)
*/

TYPE T == CHAR;

CONST
  All           == SET OF T {'\000'.. '\377'};

  Controls      == SET OF T {'\000'.. '\037', '\177'.. '\237'};
  Graphics      == All - Controls;
  /*
    The printable ISO-Latin-1 characters, including the plain 
    ASCII space '\040' and the ISO non-breaking space '\240',
    but excluding TAB, newline, etc. */

  Spaces        == SET OF T {' ', '\t', '\n', '\r', '\f', '\240'};
  /*
    Various kinds of whitespace. */

  Printing      == Graphics - Spaces;
  /*
    Characters that have a visible manifestation. */
  
  Digits        == SET OF T {'0'.. '9'};
  Uppers        == SET OF T {'A'.. 'Z', '\300'.. '\326', '\330'.. '\336'};
  Lowers        == SET OF T {'a'.. 'z', '\337'.. '\366', '\370'.. '\377'};
  Letters       == Uppers + Lowers;
  /*
    Plain and accented letters, including c-cedilla, n-tilde, 
    s-zet, ae-ligature, oe-ligature, a-ring, o-slash, and the 
    Icelandic letters. 
    
    Note that s-zet and y-umlaut have no uppercase equivalent.  Note
    also that '\327' (multiply) and '\367' (divide) are NOT included. */

  AlphaNumerics == Letters + Digits;
  
  ASCIIChars    == SET OF T {'\000'.. '\177'};

  Accented      == Letters - ASCIIChars;
  /*
    Includes c-cedilla, n-tilde, s-zet, ae-ligature, a-ring,
    o-slash, and the Icelandic letters. */

VAR /*CONST*/
  Upper:    ARRAY T OF T;
  Lower:    ARRAY T OF T;
  Control:  ARRAY T OF T;
  /* 
    These constant arrays implement character conversions (mappings):

      "Upper[c]"   == the upper-case equivalent of "c" if "c" is a symbol, else "c"
      "Lower[c]"   == the lower-case equivalent of "c" if "c" is a symbol, else "c"
      "Control[c]" == the control-shifted equivalent of "c" (i.e. "BitAnd (c, 8_037))"
                     if "c" is in Graphics , else "c"
    */

/*
  Here is the full ISO Latin-1 table:

    | 000 nul| 001 soh| 002 stx| 003 etx| 004 eot| 005 enq| 006 ack| 007 bel|
    | 010 bs | 011 ht | 012 nl | 013 vt | 014 np | 015 cr | 016 so | 017 si |
    | 020 dle| 021 dc1| 022 dc2| 023 dc3| 024 dc4| 025 nak| 026 syn| 027 etb|
    | 030 can| 031 em | 032 sub| 033 esc| 034 fs | 035 gs | 036 rs | 037 us |
    | 040 sp | 041  ! | 042  " | 043 !=| 044  $ | 045  % | 046  & | 047  ' |
    | 050  ( | 051  ) | 052  * | 053  + | 054  , | 055  - | 056  . | 057  / |
    | 060  0 | 061  1 | 062  2 | 063  3 | 064  4 | 065  5 | 066  6 | 067  7 |
    | 070  8 | 071  9 | 072  : | 073  ; | 074  < | 075  == | 076  > | 077  ? |
    | 100  @ | 101  A | 102  B | 103  C | 104  D | 105  E | 106  F | 107  G |
    | 110  H | 111  I | 112  J | 113  K | 114  L | 115  M | 116  N | 117  O |
    | 120  P | 121  Q | 122  RSym | 123  S | 124  T | 125  U | 126  V | 127  LSym |
    | 130  X | 131  Y | 132  Z | 133  [ | 134  \ | 135  ] | 136  ^ | 137  _ |
    | 140  ` | 141  a | 142  b | 143  c | 144  d | 145  e | 146  f | 147  g |
    | 150  h | 151  i | 152  j | 153  k | 154  l | 155  m | 156  n | 157  o |
    | 160  p | 161  q | 162  r | 163  s | 164  t | 165  u | 166  v | 167  w |
    | 170  x | 171  y | 172  z | 173  { | 174  | | 175  } | 176  ~ | 177 del|

    | 240 nbs| 241  ¡ | 242  ¢ | 243  £ | 244  ¤ | 245  ¥ | 246  ¦ | 247  § |
    | 250  ¨ | 251  © | 252  ª | 253  « | 254  ¬ | 255  ­ | 256  ® | 257  ¯ |
    | 260  ° | 261  ± | 262  ² | 263  ³ | 264  ´ | 265  µ | 266  ¶ | 267  · |
    | 270  ¸ | 271  ¹ | 272  º | 273  » | 274  ¼ | 275  ½ | 276  ¾ | 277  ¿ |
    | 300  À | 301  Á | 302  Â | 303  Ã | 304  Ä | 305  Å | 306  Æ | 307  Ç |
    | 310  È | 311  É | 312  Ê | 313  Ë | 314  Ì | 315  Í | 316  Î | 317  Ï |
    | 320  Ð | 321  Ñ | 322  Ò | 323  Ó | 324  Ô | 325  Õ | 326  Ö | 327  × |
    | 330  Ø | 331  Ù | 332  Ú | 333  Û | 334  Ü | 335  Ý | 336  Þ | 337  ß |
    | 340  à | 341  á | 342  â | 343  ã | 344  ä | 345  å | 346  æ | 347  ç |
    | 350  è | 351  é | 352  ê | 353  ë | 354  ì | 355  í | 356  î | 357  ï |
    | 360  ð | 361  ñ | 362  ò | 363  ó | 364  ô | 365  õ | 366  ö | 367  ÷ |
    | 370  ø | 371  ù | 372  ú | 373  û | 374  ü | 375  ý | 376  þ | 377  ÿ |

  Characters '\200' .. '\237' are unassigned.

  Some significant control codes:

    '\000' == nul == C end-of-string marker
    '\007' == bel == aubible bell
    '\010' == bs  == backspace ('\b')
    '\011' == ht  == horizontal tab ('\t')
    '\012' == nl  == newline; line separator in Unix ('\n')
    '\014' == np  == newpage; form feed ('\f')
    '\015' == cr  == carriage return ('\r')
    '\033' == esc == escape
    '\040' == sp  == normal space

    '\240' == nbs == non-breaking space
    '\255' == ­   == soft hyphen

  */

;} ISOLatin1.

/****************************************************************************/
/* (C) Copyright 1995 Universidade Estadual de Campinas (UNICAMP)           */
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
