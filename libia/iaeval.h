/* Pseudocode interpreters for standard interval arithmetic. */
/* Last edited on 2005-09-25 00:19:40 by stolfi */ 
/* Created by Jorge Stolfi 93-12-22. */

#ifndef iaeval_H
#define iaeval_H

#include <ia.h>
#include <pcode.h>

/* These procedures evaluate a pseudocode function {fcode}
  with ordinary interval arithmetic.
    
  The procedures expect the input arguments in {reg[0..N-1]}, and
  return the results in {stack[0..M-1]}, where {N}, {M} are the input
  and output arity of the function.
  
  The {reg} and {stack} arrays must be large enough to accomodate
  the execution of {fcode}. */

void ia_eval (Interval reg[], Interval stack[], pcode_instr_t fcode[]);
  /* Evaluates the function {fcode}. */

typedef struct { Interval f; Interval df; } IntervalDiff;

void ia_eval_diff (IntervalDiff reg[], IntervalDiff stack[], pcode_instr_t fcode[]);
  /* Evaluates the function {fcode}, and its derivative with respect to some 
    parameter {t}.  
    
    In any input or output quantity, the {f} field represents the
    quantity's value, and {df} represents its derivative with respect
    to {t}. */

#endif

