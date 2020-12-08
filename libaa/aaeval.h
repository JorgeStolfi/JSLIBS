/* Pseudocode interpreters for standard interval arithmetic    */
/* Created by Jorge Stolfi 93-12-22                            */

#ifndef aaeval_H
#define aaeval_H

#include <aa.h>
#include <pcode.h>

/* These procedures evaluate a pseudocode function {fcode}
  with affine arithmetic.
    
  The procedures expect the input arguments in {reg[0..N-1]}, and
  return the results in {stack[0..M-1]}, where {N}, {M} are the input
  and output arity of the function.
  
  The {reg} and {stack} arrays must be large enough to accomodate
  the execution of {fcode}. 
  
  The procedures do not clean up the AA stack; the client is responsible
  for that. (See aa_top, aa_flush, and aa_return in <aa.h>.)
  */

void aa_eval (AAP reg[], AAP stack[], pcode_instr_t fcode[]);
  /* Evaluates the function {fcode}. */

typedef struct { AAP f; AAP df; } DiffAAP;

void aa_eval_diff (DiffAAP reg[], DiffAAP stack[], pcode_instr_t fcode[]);
  /* Evaluates the function {fcode}, and its derivative with respect to some 
    parameter {t}.  
    
    In any input or output quantity, the {f} field represents the
    quantity's value, and {df} represents its derivative with respect
    to {t}. */

#endif

