/* AA test: finding the zeros of a 1-argument function. */
/* Last edited on 2005-09-25 19:51:02 by stolfi */

#ifndef fn1_zf_H
#define fn1_zf_H

#include <fn1_functions.h>

#include <aa.h>
#include <flt.h>
#include <ia.h>

#include <stdio.h>

void fn1_zf_find_and_plot_zeros(
    char *fileprefix,
    int epsformat,
    eval_fp_t eval_fp,
    eval_ia_t eval_ia,
    diff_ia_t diff_ia,
    eval_aa_t eval_aa,
    char *title,
    Interval xd,
    Interval yd,
    Float epsilon,
    Float delta,
    int m
  );
  /* 
    Generates PostScript file(s) with plots of the curve {y = f(x)},
    within the rectangle {xd} by {yd}, and with a trace of the 
    zero-finding algorithm  with relative tolerance {epsilon} and absolute 
    tolerance {delta}.
    
    If {epsformat} is 1 (true), writes two separate Encapsulated Postscript
    files (<fileprefix>-ia.ps, <fileprefix>-aa.ps), without "showpage"
    commands or captions.
    
    If {epsformat} is 0 (false), writes a single plain Postscript file
    (<fileprefix>.ps), with one graph per page, with captions.
    
    Procedures {eval_ia} and {eval_aa} should be IA and AA versions of {f};
    Procedure {diff_ia} should compute the derivative of {f} with IA.
    Parameter {m} is the number of steps to use when plotting the "true" 
    graph of {f}. 
  */

#endif







