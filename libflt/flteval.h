#ifndef flteval_H
#define flteval_H

/* Pseudocode interpreters for standard interval arithmetic */
/* Last edited on 2005-02-14 22:59:58 by stolfi  */

#include <flt.h>
#include <pcode.h>

/*
  These procedures evaluate a pseudocode function $fcode$
  with Float arithmetic, with the IEEE rounding mode flags set to 
  "nearest".
    
  The procedures expect the input arguments in $reg[0..N-1]$, and
  return the results in $stack[0..M-1]$, where $N$, $M$ are the input
  and output arity of the function.
  
  The $reg$ and $stack$ arrays must be large enough to accomodate
  the execution of $fcode$.
  
  Created by Jorge Stolfi 93-12-22.
  */

void flt_eval (Float reg[], Float stack[], pcode_instr_t fcode[]);
  /* 
    Evaluates the function $fcode$. */

typedef struct { Float f; Float df; } FloatDiff;

void flt_eval_diff (FloatDiff reg[], FloatDiff stack[], pcode_instr_t fcode[]);
  /* 
    Evaluates the function $fcode$, and its derivative with respect to some 
    parameter $t$.  
    
    In any input or output quantity, the $f$ field represents the
    quantity's value, and $df$ represents its derivative with respect
    to $t$. */

typedef struct { Float f; Float df[4]; } FloatDiff4;

void flt_eval_diff4 (FloatDiff4 reg[], FloatDiff4 stack[], pcode_instr_t fcode[]);
  /* 
    Evaluates the function $fcode$, and its partial derivatives with
    respect to four parameters $x[0..3]$.
    
    In any input or output quantity, the $f$ field represents the
    quantity's value, and $df[i]$ represents its partial derivative
    with respect to parameter $x[i]$. */

#endif

