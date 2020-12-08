#ifndef _H
#define _H


/* Miscellaneous utility functions for the "Dicio" project */
/* See the copyright and disclaimer note at the end of this file. */

IMPORT Date, Time, Rd, Wr;

PROCEDURE GetDate(): Date.T;
/* 
  Returns the current local date and time as a Date.T record. */

PROCEDURE GetUTCDate(): Date.T;
/* 
  Returns the current UTC (Greenwich) date and time as a Date.T record. */
  
PROCEDURE GetUserName(effective: BOOLEAN = TRUE): char *;
PROCEDURE GetGroupName(effective: BOOLEAN = TRUE): char *;
/* 
  Returns the user/group name of the process, or "\"\"" if not known. */

PROCEDURE GetHostName(): char *;
/* 
  Returns the name of the machine where the process is running,
  or "\"\"" if not known. */

PROCEDURE FmtDate(READONLY date: Date.T): char *;
/*
  Formats a Date.T as "Mon 1995-01-03 07:15:03 BRZ (+09:00)" */

PROCEDURE FmtDateNum(READONLY date: Date.T): char *;
/*
  Formats a Date.T as "1995-01-03-071503". */

PROCEDURE FmtTime(t: Time.T): char *;
/*
  Formats a time amount as "ddd+hh:hh:hh".  The "ddd+" is omitted
  if "t" is less than one day. */

PROCEDURE ElapsedTime(READONLY t1, t2: Date.T): Time.T;
/*
  Seconds elapsed from t1 to t2 */

/***********************************************************************/
/*  EXECUTION TIMES                                                    */
/***********************************************************************/

TYPE
  ExecTimes == RECORD
      user:   Time.T;
      system: Time.T;
      total:  Time.T;
    ;};

  ExecTimesText == RECORD
      char *user;
      char *system;
      char *total;
    ;};

PROCEDURE GetExecTimes(): ExecTimes;
/* 
  Returns user, system and total execution times since last invocation 
  of the procedure (or since the beginning of the process). */

PROCEDURE FmtExecTimes(READONLY t: ExecTimes): ExecTimesText;
/*
  Formats the times using "FmtTime". */

PROCEDURE GetExecTimesText(): ExecTimesText;
/* 
  Equivalent to "FmtExecTimes(GetExecTimes())". */

/***********************************************************************/
/*  FILE OPENING CREATION                                              */
/***********************************************************************/

PROCEDURE OpenRd(name: char *; err: Wr.T): Rd.T;
  /*
    Open file "name" for input.
    If "name" is "-", returns "Stdio.stdin". 
    If "name" is the empty string, returns NULL. 
    Prints message to "err" if the file can't be opened. */

PROCEDURE OpenWr(name: char *; err: Wr.T): Wr.T;
  /*
    Open file "name" for output. 
    If "name" is "-", returns "Stdio.stdout". 
    If "name" is the empty string, returns NULL. 
    Prints message to "err" if the file can't be opened. */

/***********************************************************************/
/*  COMMENT char *FORMATTING                                            */
/***********************************************************************/

PROCEDURE PrintDoc(wr: Wr.T; doc: char *; prefix: char *= "|");
  /*
    Writes the given "doc" text to "wr", with the "prefix" string
    in front of every line.  Supplies a final '\n' if the text is 
    non-empty but does not end with newline.*/

PROCEDURE ToLowerCase(t: char *): char *;
  /*
    Converts all uppercase ISO-Latin-1 characters in "t" to their 
    lowercase equivalents. */

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

