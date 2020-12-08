/* Functions for testing 2-dimensional zero-finders. */
/* Last edited on 2005-09-25 14:21:49 by stolfi */

#ifndef fn2_functions_H
#define fn2_functions_H

#include <aa.h>
#include <ia.h>
#include <flt.h>

/* Types of procs that evaluate a function in several ways:  */
typedef Float eval_fp_t (Float x, Float y); /* Evaluates in floating point. */
typedef Interval eval_ia_t (Interval x, Interval y); /* Evaluates with standard IA. */
typedef AAP eval_aa_t (AAP x, AAP y); /* Evaluates with affine arithmetic. */

typedef struct fn2_data_t 
  { 
    char *tag;             /* Function tag, for file names etc.. */
    char *descr;           /* Description of function. */
    eval_fp_t *eval_fp;    /* Evaluates the function in floating-point. */
    eval_ia_t *eval_ia;    /* Evaluates the function with standard IA. */
    eval_aa_t *eval_aa;    /* Evaluates the function with affine arith. */
    Interval xd, yd;       /* Interesting range of X and Y arguments. */
    
  } fn2_data_t;

fn2_data_t fn2_from_tag(char *tag);

#endif

