/* Functions for testing univariate function plotters. */
/* Last edited on 2005-09-25 14:59:57 by stolfi */

#ifndef fn1_functions_H
#define fn1_functions_H

#include <aa.h>
#include <ia.h>
#include <flt.h>

/* Types of procs that evaluate a function in several ways:  */
typedef Float eval_fp_t (Float x); /* Computes value with floating point. */
typedef Interval eval_ia_t (Interval x); /* Computes value with IA. */
typedef Interval diff_ia_t (Interval x); /* Computes derivative with IA. */
typedef AAP eval_aa_t (AAP x); /* Computes value with affine arithmetic. */

typedef struct fn1_data_t 
  { 
    char *tag;             /* Function tag, for file names etc.. */
    char *descr;           /* Description of function. */
    eval_fp_t *eval_fp;    /* Evaluates the function in floating-point. */
    eval_ia_t *eval_ia;    /* Evaluates the function with standard IA. */
    diff_ia_t *diff_ia;    /* Evaluates the derivative of the function with standard IA. */
    eval_aa_t *eval_aa;    /* Evaluates the function with affine arithmetic. */
    Interval xd;           /* Interesting range of argument. */
    Interval yd;           /* Plot range for result. */
    /* Suggested parameters for plotting or zero-finding: */
    double epsilon;        /* Relative tolerance. */
    double delta;          /* Absolute tolerance. */
    int nsub;              /* Number of sub-intervals. */
  } fn1_data_t;

fn1_data_t fn1_from_tag(char *tag);

#endif
