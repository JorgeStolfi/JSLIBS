/* Generic mathematical pseudocode */
/* Last edited on 2023-02-22 20:42:03 by stolfi */ 

#ifndef pcode_H
#define pcode_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

typedef enum {
    pcode_op_given =   0,
    pcode_op_const =   1,
    pcode_op_load =    2,
    pcode_op_store =   3,
    pcode_op_add =     4,
    pcode_op_sub =     5,
    pcode_op_neg =     6,
    pcode_op_inv =     7,
    pcode_op_mul =     8,
    pcode_op_div =     9,
    pcode_op_sqr =    10,
    pcode_op_sqrt =   11,
    pcode_op_abs =    12,
    pcode_op_max =    13,
    pcode_op_min =    14,
    pcode_op_return = 15  /* Must be the largest opcode */
  } pcode_op_t;
  
typedef struct { 
    pcode_op_t op;  /* Operation code */
    int32_t arg;        /* Operation's argument */
  } pcode_instr_t;
  
typedef struct {
    char *title;          /* Function's legible description */
    int32_t nin;          /* Number of input arguments. */
    int32_t nout;         /* Number of output results. */
    int32_t nregs;        /* Number of working registers needed. */
    int32_t nstack;       /* Number of stack entries needed. */
    int32_t nops;         /* Number of instructions */
    pcode_instr_t *code;  /* Instructions */
  } pcode_proc_t;
  
void pcode_proc_null(pcode_proc_t *pc);
  /* Sets {*pc} to a null p-code procedure, with all fields {NULL} or 0. */

pcode_proc_t pcode_parse (FILE *f);
  /*
    Parses a pseudocode procedure from $f$, up to the first RETURN
    instruction.  The input syntax is

       # first title line
       # second title line
       ...
       # last tile line
       GIVEN nin
       LOAD 0
       LOAD 2
       +
       ...
       STORE 3
       RETURN nout 
    
    The fields $nin$ and $nout$ are obtained from the "GIVEN" and
    "RETURN" directives, The fields $nops$, $nregs$, and $nstack$ are
    computed by pcode_parse itself.
    
    When a pseudocode procedure is evaluated, the arguments are passed
    as the first $nin$ registers, and the results are returned in the
    first $nout$ stack entries.
    
    At the moment, the valid opcodes are 
      
      LOAD  STORE  CONST  
      +  -  NEG  INV *  /  SQR  SQRT
      ABS MAX MIN
      
    */

void pcode_print (FILE *f, pcode_proc_t proc);
  /* 
    Prints $proc$ to the given file, according to the syntax used
    by pcode_parse. */

#endif
