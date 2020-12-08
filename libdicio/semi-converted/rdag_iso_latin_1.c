/* Last edited on 1999-06-05 20:32:58 by stolfi */

#ifndef _H
#define _H


/*
  Based on "ISOChar.m3" by Jim Meehan and Henri Gouraud 
  (Copyright 1994, Digital Equipment Corp.)
*/

IMPORT ASCII, Word;

{
  for (c = '\000' TO '\377'){ Upper [c] = ASCII.Upper [c] ;};
  for (c = '\340' TO '\376'){
    if ((c!='\367')){
      Upper [c] = VAL (ORD (c) - ORD ('a') + ORD ('A'), CHAR)
    ;}
  ;};

  for (c = '\000' TO '\377'){ Lower [c] = ASCII.Lower [c] ;};
  for (c = '\300' TO '\336'){
    if ((c!='\327')){
      Lower [c] = VAL (ORD (c) - ORD ('A') + ORD ('a'), CHAR)
    ;}
  ;};

  for (c = '\000' TO '\377'){
    if ((c IN Graphics)){
      Control [c] = VAL (Word.And (ORD (c), 8_37), CHAR)
    }else{
      Control [c] = c
    ;}
  ;}
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
