#ifndef _H
#define _H


/* See the copyright and disclaimer note at the end of this file. */

IMPORT Rd, Wr, TextF, TextWr, Fmt, Thread; 
IMPORT Basics;

FROM Basics IMPORT Symbol, String, Done;
FROM Basics IMPORT NAT, CHARS;
FROM Encoding IMPORT BadChar, BadChars, BadLetter, BadString;

REVEAL 
  T == Public BRANDED 
    OBJECT 
    OVERRIDES
      StringToChars = StringToChars;
      StringToText = StringToText;
      PrintString = PrintString;
      PrintLetter = PrintLetter;

      CharsToString = CharsToString;
      TextToString = TextToString;
      ReadString = ReadString;

      CharToLetter = CharToLetter;
      SymbolToChar = SymbolToChar;
    ;};

PROCEDURE StringToChars(
    T e;
    String READONLY w; 
    VAR /*IO*/ c: REF CHARS; 
    VAR /*OUT*/ len: NAT;
  ) RAISES {BadString} ==
  {
    len = NUMBER(w);
    if ((c == NULL) || (NUMBER(c^) < len)){ Basics.ExpandChars(c, len) ;};
    with (cc == c^){
      for (i = 0 TO LAST(w)){ 
        TRY
          cc[i] = SymbolToChar(e, w[i])
        EXCEPT
        | BadLetter ==> RaiseBadString(w[i], i)
        ;} 
      ;}
    ;}
  ;} StringToChars;

PROCEDURE StringToText(
    T e;
    String READONLY w;
  ): char *RAISES {BadString} ==
  {
    with (wr == TextWr.New()){
      PrintString(e, wr, w);
      return TextWr.ToText(wr)
    ;}
  ;} StringToText;

PROCEDURE PrintString(
    <*UNUSED); e: T;
    wr: Wr.T;
    String READONLY w;
  ) RAISES {BadString} ==
  <* FATAL Wr.Failure, Thread.Alerted );
  {
    for (i = 0 TO LAST(w)){ 
      with (wi == w[i]){
        if ((NOT (wi IN ValidLetters))){ RaiseBadString(wi, i) ;};
        Wr.PutChar(wr, VAL(wi, CHAR))
      ;}
    ;}
  ;} PrintString;

PROCEDURE RaiseBadString(let: Symbol; i: NAT) RAISES {BadString} ==
  {
    RAISE BadString(
      "bad symbol at position " & Fmt.Int(i) &
      " == " & Fmt.Int(let) & " == 8_" & Fmt.Int(let, base = 8)
    )
  ;} RaiseBadString;

CONST OctalDigit == ARRAY [0..7] OF CHAR{'0', '1', '2', '3', '4', '5', '6', '7'};

PROCEDURE PrintLetter(
    T e;
    wr: Wr.T; 
    symbol: Symbol;
  ) ==
  <* FATAL Wr.Failure, Thread.Alerted );
  {
    TRY
      Wr.PutChar(wr, SymbolToChar(e, symbol))
    EXCEPT
    | BadLetter(bad) ==>
        Wr.PutChar(wr, '\\');
        with (n == ORD(bad)){
          Wr.PutChar(wr, OctalDigit[n DIV 64]);          
          Wr.PutChar(wr, OctalDigit[(n DIV 8) MOD 8]);          
          Wr.PutChar(wr, OctalDigit[n MOD 8]);          
        ;};
    ;};
  ;} PrintLetter;

PROCEDURE CharsToString(
    T e;
    CHARS READONLY c;
    VAR /*IO*/ w: REF String;
    VAR /*OUT*/ len: NAT;
  ) RAISES {BadChars} ==
  {
    len = NUMBER(c);
    if ((w == NULL) || (NUMBER(w^) < len)){ Basics.ExpandString(w, len) ;};
    with (ww == w^){
      for (i = 0 TO LAST(c)){ 
        TRY
          ww[i] = CharToLetter(e, c[i])
        EXCEPT
        | BadChar ==> RaiseBadChars(c[i], i)
        ;} 
      ;}
    ;}
  ;} CharsToString;

PROCEDURE TextToString(
    T e;
    char *t;
    VAR /*IO*/ w: REF String;
    VAR /*OUT*/ len: NAT;
  ) RAISES {BadChars} ==
  {
    with (tt == t^){
      assert(tt[LAST(tt)] == '\000' );
      CharsToString(e, SUBARRAY(tt, 0, NUMBER(tt)-1), w, len)
    ;};
  ;} TextToString;
  
CONST SpaceLetter: Symbol == ORD(' ');

PROCEDURE ReadString(
    T e;
    rd: Rd.T;
    VAR /*IO*/ w: REF String;
    VAR /*OUT*/ len: NAT;
  ) RAISES {BadChars, Done, Rd.Failure, Thread.Alerted} ==
  CHAR *c;
  {
    len = 0;
    TRY 
      REPEAT c = Rd.GetChar(rd) UNTIL c!=' '
    EXCEPT 
      Rd.EndOfFile ==> RAISE Done 
    ;};
    while (c!='\n'){
      Basics.ExpandString(w, len+1);
      TRY
        w[len] = CharToLetter(e, c)
      EXCEPT 
      | BadChar ==> RaiseBadChars(c, len)
      ;};
      INC(len);
      TRY 
        c = Rd.GetChar(rd)
      EXCEPT 
        Rd.EndOfFile ==> c = '\n'
      ;};
    ;};
    while (len > 0)  AND  AND  (w[len-1] == SpaceLetter){ DEC(len) ;};
  ;} ReadString;

PROCEDURE RaiseBadChars(c: CHAR; i: NAT) RAISES {BadChars} ==
  {
    RAISE BadChars(
      "bad character at position " & Fmt.Int(i) &
      " == \\" & Fmt.Int(ORD(c), base = 8)
    )
  ;} RaiseBadChars;

PROCEDURE CharToLetter(
    <*UNUSED); e: T; 
    CHAR ch;
  ): Symbol RAISES {BadChar} ==
  {
    if ((NOT(ch IN ValidChars))){ RAISE BadChar(ch) ;};
    return ORD(ch)
  ;} CharToLetter;

PROCEDURE SymbolToChar(
    <*UNUSED); e: T; 
    let: Symbol;
  ): CHAR RAISES {BadLetter} ==
  {
    if ((NOT(let IN ValidLetters))){ RAISE BadLetter(let) ;};
    return VAL(let, CHAR)
  ;} SymbolToChar;

PROCEDURE Init(<*UNUSED); e: T) ==
  {
  ;} Init;

PROCEDURE New(): T ==
  {
    with (e == NEW(T)){ 
      Init(e);
      return e
    ;};
  ;} New;

{
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
