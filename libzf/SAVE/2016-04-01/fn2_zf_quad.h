/* Plots of zeros of a 2-ary function with adaptive k-d-tree. */
/* Last edited on 2005-09-25 21:29:26 by stolfi */ 

#ifndef fn2_zf_quad_H
#define fn2_zf_quad_H

#include <fn2_functions.h>

#include <aa.h>
#include <ia.h>
#include <flt.h>

#include <stdio.h>

void fn2_zf_quad_plots(
    char *fileprefix,
    int epsformat,
    char *title,
    eval_fp_t *eval_fp,    /* Evaluates the function in floating-point. */
    eval_ia_t *eval_ia,    /* Evaluates the function with standard IA. */
    eval_aa_t *eval_aa,    /* Evaluates the function with affine arith. */
    Interval xd,
    Interval yd,
    int n,
    int m
  );
  /* 
    Generates PostScript file(s) with plots of a curve {f(x,y)=0}, 
    within the rectangle {xd} by {yd}.
  
    If {epsformat} is 1 (true), writes two separate Encapsulated Postscript
    files (<fileprefix>-ia.ps, <fileprefix>-aa.ps), without "showpage"
    commands or captions.
    
    If {epsformat} is 0 (false), writes a single plain Postscript file
    (<fileprefix>.ps), with one graph per page, with captions.
    
    {eval_ia} and {eval_aa} should be Interval and AA versions of {f}.
    {n} is the grid size used for testing these two versions.
    {m} is the grid size used for plotting the "true" zeros of {f}. */

#endif




