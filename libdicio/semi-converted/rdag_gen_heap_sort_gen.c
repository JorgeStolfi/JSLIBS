#ifndef _H
#define _H


/* See the copyright and disclaimer note at the end of this file. */

PROCEDURE Sort(
    VAR /*IO*/ s: ARRAY OF X.T; 
    VAR /*IO*/ ss: REF ARRAY OF Y.T;
  ) == 
  unsigned *q, r, mx;
       cc,smx: X.T;
       ccss: Y.T;
  {

    with (
      n == NUMBER(s),
      tss == ss#NULL
   ){
      /* 1. Build heap with largest element at s[0] */
      for ( p = 0  TO n-1 ){
        cc = s[p]; 
        if ((tss)){ ccss = ss[p] ;};
        r = p;
        while (1){
          if (( r==0 )){  EXIT  ;};
          q = (r-1) DIV 2;
          if ((X.Compare(s[q],cc) == -1 )){
            s[r] = s[q]; 
            if ((tss)){ ss[r] = ss[q] ;};
          }else{
            EXIT
          ;};
          r = q
        ;};
        s[r] = cc; 
        if ((tss)){ ss[r] = ccss ;}
      ;};
      /* 2. Remove elements from heap and insert at end */
      for ( p = n-1  TO 1  BY -1 ){
        /* save s[p] */
        cc = s[p]; 
        if ((tss)){ ccss = ss[p] ;};
        /* Move largest heap element to pos[p] */
        s[p] = s[0]; 
        if ((tss)){ ss[p] = ss[0] ;};
        /* Insert cc in remaining heap, from root down  */
        q = 0;
        while (1){
          s[q] = cc; 
          if ((tss)){ ss[q] = ccss ;};
          r = 2*q+1;
          /* Find largest among cc, s[LEFT(q)], s[RIGHT(q)]  */
          mx = q; smx = cc;
          if (((r<p))  AND  AND  ((X.Compare(s[r],smx) == +1) )){
            mx = r; smx = s[r]
          ;};
          INC(r);
          if (( (r<p))  AND  AND  ((X.Compare(s[r],smx) == +1) )){
            mx = r; smx = s[r]
          ;};
          /*  See who won  */
          if (( mx==q )){
            /*  Stop here  */
            EXIT
          }else{
            /*  Promote child and advance  */
            s[q] = s[mx];  
            if ((tss)){ ss[q] = ss[mx] ;}; 
            q = mx
          ;}
        ;}
      ;}
    ;}
  ;} Sort;

{
;} GenHeapSort.

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
 
