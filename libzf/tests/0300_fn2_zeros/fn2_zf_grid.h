/* Plots of zeros of a 2-ary function with uniform grid. */
/* Last edited on 2023-02-18 09:24:21 by stolfi */ 

#ifndef fn2_zf_grid_H
#define fn2_zf_grid_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <aa.h>
#include <flt.h>
#include <ia.h>

#include <fn2_functions.h>

void fn2_zf_grid_plots(
    char *prefix,
    char *title,
    eval_fp_t *eval_fp, /* Evaluates the function in floating-point. */
    eval_ia_t *eval_ia, /* Evaluates the function with standard IA. */ 
    eval_aa_t *eval_aa, /* Evaluates the function with affine arith. */
    Interval xd,
    Interval yd,
    int32_t n,
    int32_t m
  );
  /* Generates plots of the implicit curve {f(x,y)=0}, within the
   rectangle {xd} by {yd}, by evaluating {f} with IA or AA on cells of a
   uniform grid with {n} rows and {n} columns.
    
    Writes two EPS files named "out/{prefix}_ia.eps" and
    "out/{prefix}_aa.eps" 
    
    The procedure {eval_fp} should evaluate {f(x,y)} given two {Float}
    arguments {x,y}, with floating point. It is used to draw the "true"
    plot of {f(x,y)=0}, obtained by evaluating {f} on a finer grid of
    {m} rows and columns.
    
    The procedures {eval_ia} and {eval_aa} must compute {f(x,y)} given
    intervals {xd,yd} for {x} and {y}, with IA (returning an {Interval})
    and with AA (returning an affine form). */

#endif







