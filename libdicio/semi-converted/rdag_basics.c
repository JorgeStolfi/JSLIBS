#ifndef _H
#define _H


/* See the copyright and disclaimer note at the end of this file. */

PROCEDURE ExpandString(VAR s: REF String; n: unsigned) ==
  {
    with (nold == NUMBER(s^)){
      if ((nold < n)){
        with (
          r == NEW(REF String, MAX(n + 10, 2*nold))
       ){
          SUBARRAY(r^, 0, nold) = s^;
          s = r
        ;}
      ;}
    ;}
  ;} ExpandString;

PROCEDURE ExpandChars(VAR s: REF CHARS; n: unsigned) ==
  {
    with (nold == NUMBER(s^)){
      if ((nold < n)){
        with (
          r == NEW(REF CHARS, MAX(n + 10, 2*nold))
       ){
          SUBARRAY(r^, 0, nold) = s^;
          s = r
        ;}
      ;}
    ;}
  ;} ExpandChars;

PROCEDURE CopyString(READONLY s: String): REF String ==
  {
    with (r == NEW(REF String, NUMBER(s))){
      r^ = s;
      return r
    ;}
  ;} CopyString;

{
;} Basics.

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

