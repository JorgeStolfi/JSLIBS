/* Plotting graph of a 1-argument function with affine arithmetic. */
/* Last edited on 2005-09-25 18:55:49 by stolfi */

#ifndef allgraphs_H
#define allgraphs_H

#include <fn1_functions.h>

#include <flt.h>
#include <aa.h>
#include <ia.h>

#include <stdio.h>

void allgraphs_plot(
    char *fileprefix,
    int epsformat,
    char *title,
    eval_fp_t eval_fp,
    eval_ia_t eval_ia,
    diff_ia_t diff_ia,
    eval_aa_t eval_aa,
    Float xmin, Float xmax,
    Float ymin, Float ymax,
    int nsub,
    int nsteps
  );
  /* 
    Generates PostScript file(s) with plots of the graph of a function
    {F(x)}, on the interval {[xmin .. xmax]} computed with floating point,
    ordinary interval arithmetic and with affine arithmetic.
    
    If {epsformat} is 1 (true), writes two separate Encapsulated Postscript
    files (<fileprefix>-ia.ps, <fileprefix>-aa.ps), without "showpage"
    commands or captions.
    
    If {epsformat} is 0 (false), writes a single plain Postscript file
    (<fileprefix>.ps), with one graph per page, with captions.
    
    The functions {eval_ia} and {eval_aa} should be IA and AA versions of {f}.
    {xd} and {yd} are the plot coordinate ranges. {nsub} is the number of
    intervals into which {xd} is to be divided. {nsteps} is the number of
    steps for the floating-point graph of {F}. */

#endif
