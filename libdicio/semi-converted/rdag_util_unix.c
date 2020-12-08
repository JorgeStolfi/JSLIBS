/* Last edited on 1999-06-05 20:50:38 by stolfi */ 

#ifndef _H
#define _H


/* See the copyright and disclaimer note at the end of this file. */

IMPORT M3toC, Unix, Uugid, Upwd, Ugrp, Fmt;

PROCEDURE GetUserName(effective: BOOLEAN = TRUE): char *==
  INTEGER *uid;
  {
    if ((effective)){ 
      uid = Uugid.geteuid()
    }else{
      uid = Uugid.getuid()
    ;};
    with (user == (Upwd.getpwuid(uid)).pw_name){
      if ((user == NULL)){
        return Fmt.Int(uid)
      }else{
        return M3toC.CopyStoT(user)
      ;}
    ;}
  ;} GetUserName;

PROCEDURE GetGroupName(effective: BOOLEAN = TRUE): char *==
  INTEGER *gid;
  {
    if ((effective)){ 
      gid = Uugid.getegid()
    }else{
      gid = Uugid.getgid()
    ;};
    with (group == (Ugrp.getgrgid(gid)).gr_name){
      if ((group == NULL)){
        return Fmt.Int(gid)
      }else{
        return M3toC.CopyStoT(group)
      ;}
    ;}
  ;} GetGroupName;

PROCEDURE GetHostName(): char *== 
  CONST MaxHostNameChars == 256;
  VAR buf: ARRAY [0..MaxHostNameChars+1] OF CHAR;
  {
    with (res == Unix.gethostname(ADR(buf[0]), NUMBER(buf)-1)){
      if ((res!=0)){ return "" ;};
      buf[LAST(buf)] = '\000';
      return M3toC.CopyStoT(buf)
    ;};
  ;} GetHostName;

{
;} UtilUnix.

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
