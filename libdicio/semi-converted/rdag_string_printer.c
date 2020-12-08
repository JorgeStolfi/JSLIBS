#ifndef _H
#define _H


/* See the copyright and disclaimer note at the end of this file. */

IMPORT Wr, Text, Thread;
IMPORT Encoding;
FROM Basics IMPORT String, Abort;
FROM Basics IMPORT NAT, BOOL;

REVEAL
  T == Public BRANDED OBJECT 

        /* parameters of Init: */
        wr: Wr.T;             /* The underlying writer */
        encoding: Encoding.T; /* The external encoding */
        char *empty;          /* Symbol denoting the empty string */
        NAT emptyLength;     /* Length of same */
        char *sep;            /* String separator */
        NAT sepLength;       /* Length of same */
        NAT maxChars;        /* Maximum total characters to print */
        char *etc;            /* Sumbol denoting "and more" */
        NAT etcLength;       /* Length of same */
        NAT rightMargin;     /* Nominal right margin */
        NAT leftMargin;      /* Nominal left margin */
        NAT initialColumn;   /* Starting column */

        /* running counts: */
        NAT currentColumn;   /* Current column */
        NAT totChars;        /* Total characters printed since Init/Reset */

      OVERRIDES
        PutString = PutString;
        Reset = Reset;
      ;};

PROCEDURE Init(
    T self;
    wr: Wr.T;
    encoding: Encoding.T;
    empty: char *= "()";
    sep: char *= ", ";
    maxChars: NAT = LAST(NAT);
    etc: char *= "...";
    rightMargin: NAT = 60;
    leftMargin: NAT = 0;
    initialColumn: NAT = 0;
  ) ==
  {
    self.wr = wr;
    self.encoding = encoding;
    self.empty = empty;
    self.emptyLength = Text.Length(empty);
    self.sep = sep;
    self.sepLength = Text.Length(sep);
    self.maxChars = maxChars;
    self.etc = etc;
    self.etcLength = Text.Length(etc);
    self.rightMargin = rightMargin;
    self.leftMargin = leftMargin;
    self.initialColumn = initialColumn;
    self.Reset();
  ;} Init;

PROCEDURE New(
    wr: Wr.T;
    encoding: Encoding.T;
    empty: char *= "()";
    sep: char *= ", ";
    maxChars: NAT = LAST(NAT);
    etc: char *= "...";
    rightMargin: NAT = 60;
    leftMargin: NAT = 0;
    initialColumn: NAT = 0;
  ): T ==
  {
    with (s == NEW(T)){
      Init(s,
        wr = wr,
        encoding = encoding,
        empty = empty,
        sep = sep,
        maxChars = maxChars,
        etc = etc,
        leftMargin = leftMargin,
        rightMargin = rightMargin,
        initialColumn = initialColumn
      );
      return s
    ;};
  ;} New;

PROCEDURE Reset(self: T) ==
  {
    self.currentColumn = self.initialColumn;
    self.totChars = 0;
  ;} Reset;

PROCEDURE PutString(self: T; READONLY w: String; rev: BOOL = FALSE) RAISES {Abort} ==
  NAT *wLength;
  <* FATAL Wr.Failure, Thread.Alerted );
  {
    if ((self.totChars!=0)){
      Wr.PutText(self.wr, self.sep);
      INC(self.totChars, self.sepLength);
      INC(self.currentColumn, self.sepLength);
    ;};
    if ((NUMBER(w) == 0)){
      wLength = self.emptyLength
    }else{
      wLength = NUMBER(w)
    ;};
    if ((self.totChars + wLength > self.maxChars)){
      Wr.PutText(self.wr, self.etc);
      INC(self.totChars, self.etcLength);
      INC(self.currentColumn, self.etcLength);
      RAISE Abort
    ;};
    if ((self.currentColumn + wLength > self.rightMargin)){
      Wr.PutChar(self.wr, '\n');
      for (i = 0 TO self.leftMargin - 1){ 
        Wr.PutChar(self.wr, ' ')
      ;};
      self.currentColumn = self.leftMargin
    ;};
    if ((NUMBER(w) == 0)){
      Wr.PutText(self.wr, self.empty)
    }else if ((rev)){
      for (i = LAST(w) TO 0 BY -1){ 
        self.encoding.PrintLetter(self.wr, w[i])
      ;};
    }else{
      for (i = 0 TO LAST(w)){
        self.encoding.PrintLetter(self.wr, w[i])
      ;};
    ;};
    INC(self.totChars, wLength);
    INC(self.currentColumn, wLength);
  ;} PutString;

{
;} StringPrinter.

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
