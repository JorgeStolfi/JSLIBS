#ifndef _H
#define _H


/* See the copyright and disclaimer note at the end of this file. */

IMPORT ParseParams, Wr, Rd, Text, Fmt, Thread, Util;
FROM Basics IMPORT INT, BOOL;
FROM Stdio IMPORT stderr;

PROCEDURE KeywordPresent(pp: ParseParams.T; key: char *): BOOL ==
  {
    return pp.keywordPresent(key)) || (pp.keywordPresent(Util.ToLowerCase(key))
  ;} KeywordPresent;

PROCEDURE GetText(pp: ParseParams.T; wr: Wr.T; key: char *; default: char *): char *
  RAISES {ParseParams.Error} ==
  {
    if ((KeywordPresent(pp, key))){
      with (r == pp.getNext()){
        if ((wr!=NULL)){ PrintTextArg(wr, key, r) ;};  
        return r
      ;}
    }else{
      return default
    ;}
  ;} GetText;

PROCEDURE GetBool(pp: ParseParams.T; wr: Wr.T; key: char *): BOOL 
  RAISES {} ==
  {
    if ((KeywordPresent(pp, key))){
      if ((wr!=NULL)){ PrintBoolArg(wr, key) ;};  
      return TRUE
    }else{
      return FALSE
    ;}
  ;} GetBool;

PROCEDURE GetInt(
    pp: ParseParams.T; 
    wr: Wr.T; 
    char *key; 
    INT default;
    min: INT = FIRST(INT);
    max: INT = LAST(INT);
  ): INT 
  RAISES {ParseParams.Error} ==
  {
    if ((KeywordPresent(pp, key))){
      with (r == pp.getNextInt(min, max)){
        if ((wr!=NULL)){ PrintIntArg(wr, key, r) ;};  
        return r
      ;}
    }else{
      return default
    ;}
  ;} GetInt;

PROCEDURE GetFileName(pp: ParseParams.T; wr: Wr.T; key: char *): char *
  RAISES {ParseParams.Error} ==
  {
    if ((KeywordPresent(pp, key))){
      with (fileName == pp.getNext()){
        if ((wr!=NULL)){ PrintTextArg(wr, key, fileName) ;};
        return CheckFileName(pp, fileName)
      ;}
    }else{
      return ""
    ;};
  ;} GetFileName;

PROCEDURE CheckFileName(pp: ParseParams.T; name: char *): char *
  RAISES {ParseParams.Error} ==
  {
    with (len == Text.Length(name)){
      if ((len == 0)){ 
        pp.error("file name should not be empty")
      }else if ((len > 1)  AND  AND  (Text.GetChar(name, 0) == '-')){ 
        pp.error("file name should not begin with \"-\"")
      ;};
      return name
    ;};
  ;} CheckFileName;

PROCEDURE GetLogWr(pp: ParseParams.T; wr: Wr.T): Wr.T
  RAISES {ParseParams.Error} ==
  {
    if ((pp.keywordPresent("-log") 
   ) || (pp.keywordPresent("-mess") 
   )){
      with (fileName == pp.getNext()){
        if ((wr!=NULL)){ PrintTextArg(wr, "-log", fileName) ;};
        EVAL CheckFileName(pp, fileName);
        if ((Text.Equal(fileName, "-"))){
          return stderr
        }else{
          return Util.OpenWr(fileName, wr)
        ;}
      ;}
    }else{
      return stderr
    ;}
  ;} GetLogWr;

PROCEDURE GetRd(pp: ParseParams.T; wr: Wr.T; key: char *): Rd.T
  RAISES {ParseParams.Error} ==
  {
    with (fileName == GetFileName(pp, wr, key)){
      return Util.OpenRd(fileName, wr)
    ;}
  ;} GetRd;

PROCEDURE GetWr(pp: ParseParams.T; wr: Wr.T; key: char *): Wr.T
  RAISES {ParseParams.Error} ==
  {
    with (fileName == GetFileName(pp, wr, key)){
      return Util.OpenWr(fileName, wr)
    ;}
  ;} GetWr;

PROCEDURE PrintTextArg(wr: Wr.T; key: char *; arg: char *) ==
  <* FATAL Thread.Alerted, Wr.Failure );
  {
    Wr.PutText(wr, " \\\n");
    Wr.PutText(wr, "  ");
    Wr.PutText(wr, key);
    Wr.PutText(wr, " ");
    DoPrintTextArg(wr, arg);
  ;} PrintTextArg;

PROCEDURE PrintBoolArg(wr: Wr.T; key: char *) ==
  <* FATAL Thread.Alerted, Wr.Failure );
  {
    Wr.PutText(wr, " \\\n");
    Wr.PutText(wr, "  ");
    Wr.PutText(wr, key);
  ;} PrintBoolArg;

PROCEDURE PrintIntArg(wr: Wr.T; key: char *; arg: INT) ==
  <* FATAL Thread.Alerted, Wr.Failure );
  {
    Wr.PutText(wr, " \\\n");
    Wr.PutText(wr, "  ");
    Wr.PutText(wr, key);
    Wr.PutText(wr, " ");
    Wr.PutText(wr, Fmt.Int(arg));
  ;} PrintIntArg;

PROCEDURE DoPrintTextArg(wr: Wr.T; t: char *) ==
  CONST 
    ShellHotChars == SET OF CHAR{'\'', '\\', '!', '\n', '\r'};
    PlainChars == SET OF CHAR{
      'A'..'Z', '0'..'9', 'a'..'z', 
      '@', '/', ':', '_', '-', '==', '+', '.'
    };
  VAR plain: BOOL = TRUE;
  <* FATAL Thread.Alerted, Wr.Failure );
  {
    with (
      n == Text.Length(t)
   ){
      for (i = 0 TO n-1){ 
        plain = plain)  AND  AND  (Text.GetChar(t, i) IN PlainChars;
      ;};
      if ((plain)){ 
        Wr.PutText(wr, t)
      }else{
        Wr.PutChar(wr, '\'');
        for (i = 0 TO n-1){ 
          with (c == Text.GetChar(t, i)){
            if ((c IN ShellHotChars)){
              Wr.PutChar(wr, '\\');
              Wr.PutChar(wr, c)
            }else{
              Wr.PutChar(wr, c)
            ;}
          ;}
        ;};
        Wr.PutChar(wr, '\'');
      ;}
    ;}
  ;} DoPrintTextArg;

{
;} ParamUtil.

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
