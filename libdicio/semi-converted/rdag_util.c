#ifndef _H
#define _H


/* See the copyright and disclaimer note at the end of this file. */

IMPORT Date, Time, Text, Fmt, Utime, ISOLatin1;
IMPORT Rd, FileRd, TextRd, Wr, FileWr, TextWr, OSError, Process, Thread;
FROM Stdio IMPORT stdin, stdout;

PROCEDURE GetDate(): Date.T ==
  {
    return Date.FromTime(Time.Now())
  ;} GetDate;
  
PROCEDURE GetUTCDate(): Date.T ==
  {
    return Date.FromTime(Time.Now(), z = Date.UTC)
  ;} GetUTCDate;

PROCEDURE ElapsedTime(READONLY d1, d2: Date.T): Time.T ==
  <* FATAL Date.Error );
  {
    return Date.ToTime(d2) - Date.ToTime(d1)
  ;} ElapsedTime;
  
PROCEDURE FmtDate(READONLY d: Date.T): char *==
  CONST
    WeekDayName == ARRAY Date.WeekDay OF char *{
      "Sun","Mon","Tue","Wed","Thu", "Fri","Sat"
    };
  {
    with (
      month == ORD(d.month)+1
   ){
      return 
        WeekDayName[d.weekDay] & " " & 
        Fmt.Pad(Fmt.Int(d.year),   4) & "-" &
        Fmt.Pad(Fmt.Int(month),    2, '0') & "-" &
        Fmt.Pad(Fmt.Int(d.day),    2, '0') & " " &
        Fmt.Pad(Fmt.Int(d.hour),   2, '0') & ":" &
        Fmt.Pad(Fmt.Int(d.minute), 2, '0') & ":" &
        Fmt.Pad(Fmt.Int(d.second), 2, '0') & " " &
        d.zone & " (" &
        FmtZoneOffset(d.offset) & ")"
    ;}
  ;} FmtDate;
  
PROCEDURE FmtZoneOffset(secs: INTEGER): char *==
/*
  Formats a timezone offset (given in seconds) as "shh:mm",
  where "s" is the sign. */
  char **sign;
  {
    if ((secs < 0)){ sign = "-" }else{ sign = "+" ;};
    return 
      sign &
      Fmt.Pad(Fmt.Int(secs DIV 3600), 2, '0') & ":" &
      Fmt.Pad(Fmt.Int((secs DIV 60) MOD 60), 2, '0')
  ;} FmtZoneOffset;
  
PROCEDURE FmtDateNum(READONLY d: Date.T): char *==
  {
    with (
      month == ORD(d.month)+1
   ){
      return 
        Fmt.Pad(Fmt.Int(d.year),        4, '0') & "-" &
        Fmt.Pad(Fmt.Int(month),         2, '0') & "-" &
        Fmt.Pad(Fmt.Int(d.day),         2, '0') & "-" &
        Fmt.Pad(Fmt.Int(d.hour),        2, '0') &
        Fmt.Pad(Fmt.Int(d.minute),      2, '0') &
        Fmt.Pad(Fmt.Int(d.second),      2, '0')
    ;}
  ;} FmtDateNum;
  
PROCEDURE FmtTime(t: Time.T): char *==
  char **d;
  {
    with (
      secs == ROUND(t),
      mins == secs DIV 60,
      hurs == mins DIV 60,
      days == hurs DIV 24
   ){
      if ((days == 0)){
        d = ""
      }else{
        d = Fmt.Pad(Fmt.Int(days), 3) & "+"
      ;};
      return 
        d &
        Fmt.Pad(Fmt.Int(hurs MOD 24), 2, '0') & ":" &
        Fmt.Pad(Fmt.Int(mins MOD 60), 2, '0') & ":" &
        Fmt.Pad(Fmt.Int(secs MOD 60), 2, '0');
    ;}
  ;} FmtTime;

VAR
  StartingTimes = NEW(Utime.struct_tms_star);
  CurrentTimes  = NEW(Utime.struct_tms_star);

PROCEDURE GetExecTimes(): ExecTimes ==
  ExecTimes *res;
  {
    with (retCode == Utime.times(CurrentTimes)){
      assert(retCode == 0 );
    ;};
    with (
      cu == CurrentTimes^.tms_utime,
      su == StartingTimes^.tms_utime,
      cs == CurrentTimes^.tms_stime,
      ss == StartingTimes^.tms_stime,
      u == cu - su,
      s == cs - ss
   ){
      res.user = ((double)u);
      res.system = ((double)s);
      res.total = ((double)u+s);
      su = cu;
      ss = cs;
      return res
    ;}
  ;} GetExecTimes;

PROCEDURE FmtExecTimes(READONLY t: ExecTimes): ExecTimesText ==
  {
    return ExecTimesText{
      user   = FmtTime(t.user),
      system = FmtTime(t.system),
      total  = FmtTime(t.total)
    }
  ;} FmtExecTimes;

PROCEDURE GetExecTimesText(): ExecTimesText ==
  {
    return FmtExecTimes(GetExecTimes())
  ;} GetExecTimesText;
  
PROCEDURE OpenRd(name: char *; err: Wr.T): Rd.T ==
  <* FATAL Thread.Alerted, Wr.Failure );
  {
    if ((Text.Empty(name))){
      return NULL
    }else if ((Text.Equal(name, "-"))){
      return stdin
    }else{
      TRY
        with (rd == FileRd.Open(name)){
          return rd
        ;}
      EXCEPT
      | OSError.E ==>
          Wr.PutText(err, "*** Util.OpenRd: error opening file == \"");
          Wr.PutText(err, name);
          Wr.PutText(err, "\"\n");
          Process.Exit(1);
          assert(FALSE );
      ;}
    ;}
  ;} OpenRd;

PROCEDURE OpenWr(name: char *; err: Wr.T): Wr.T ==
  <* FATAL Thread.Alerted, Wr.Failure );
  {
    if ((Text.Empty(name))){
      return NULL
    }else if ((Text.Equal(name, "-"))){
      return stdout
    }else{
      TRY
        with (wr == FileWr.Open(name)){
          return wr
        ;}
      EXCEPT
      | OSError.E ==>
          Wr.PutText(err, "*** Util.OpenWr: error opening file == \"");
          Wr.PutText(err, name);
          Wr.PutText(err, "\"\n");
          Process.Exit(1);
          assert(FALSE );
      ;}
    ;}
  ;} OpenWr;

PROCEDURE PrintDoc(wr: Wr.T; doc: char *; prefix: char *= "|") ==

  VAR rd: Rd.T = TextRd.New(doc);
  
  PROCEDURE CopyLine() RAISES {Rd.EndOfFile} ==
  /*
    Copy one line from "rd" to "wr", prefixed by "prefix". 
    Supplies a final '\n' if next line exists but does not end with newline.
    Raises Rd.EndOfFile if there are no more lines in "rd". */
    
    <* FATAL Rd.Failure, Wr.Failure, Thread.Alerted );
    CHAR *c;
    {
      c = Rd.GetChar(rd); /* If EOF here, propagate to caller */
      Wr.PutText(wr, prefix);
      Wr.PutChar(wr, c);
      while (c!='\n'){
        TRY c = Rd.GetChar(rd) EXCEPT Rd.EndOfFile ==> c = '\n' ;};
        Wr.PutChar(wr, c)
      ;}
    ;} CopyLine;

  {
    TRY while (1){CopyLine() ;} EXCEPT Rd.EndOfFile ==> /* Ok */ ;};
  ;} PrintDoc;

PROCEDURE ToLowerCase(t: char *): char *==
  <* FATAL Rd.EndOfFile, Rd.Failure, Wr.Failure, Thread.Alerted );
  {
    with (wr == NEW(TextWr.T).init(), rd == NEW(TextRd.T).init(t)){
      while (NOT Rd.EOF(rd)){ Wr.PutChar(wr, ISOLatin1.Lower[Rd.GetChar(rd)]) ;};
      return TextWr.ToText(wr)
    ;}
  ;} ToLowerCase;

{
  with (retCode == Utime.times(StartingTimes)){
    assert(retCode == 0 );
  ;}
;} Util.

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
