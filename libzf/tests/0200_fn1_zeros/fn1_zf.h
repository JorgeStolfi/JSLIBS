/* AA test: finding the zeros of a 1-argument function. */
/* Last edited on 2024-12-05 10:41:10 by stolfi */

#ifndef fn1_zf_H
#define fn1_zf_H

#include <stdint.h>
#include <stdio.h>

#include <ia.h>
#include <flt.h>
#include <aa.h>

#include <fn1_functions.h>

void fn1_zf_find_and_plot_zeros(
    char *prefix,
    eval_fp_t eval_fp,
    eval_ia_t eval_ia,
    diff_ia_t diff_ia,
    eval_aa_t eval_aa,
    char *title,
    Interval xd,
    Interval yd,
    double epsilon,
    double delta,
    int32_t m
  );
  /* 
    Generates Encapsulated PostScript files with plots of a curve {y =
    f(x)}, within the rectangle {xd} by {yd}, and with a trace of the
    zero-finding algorithm using different validated artithmetic
    models.
    
    The domain interval is adaptively and recursively divided into sub-intervals
    and a validated arithmetic model is used to determine the sign of
    the function in each sub-interval.  The process ends when the 
    zeros have been identified with relative tolerance {epsilon} and
    absolute tolerance {delta}.
    
    The procedure {eval_fp} should evaluate {f} at a given {Float}
    argument with floating point. Procedures {eval_ia} and {eval_aa}
    should compute {f} over a given sub-interval with {IA} and {AA},
    respectively; Procedure {diff_ia} should compute the derivative of
    {f} over an interval, with IA.
    
    The procedure writes separate Postscript files
    "out/{prefix}_{arith}_{suffix}.eps" where "{arith}" specifies the
    validated arithmetic model used:
    
      "ia" Box enclosures obtained by computing the function {f}
           with classical Interval Arithmetic.
      
      "ar" Box enclosured obtained by computing {f}
           with Affine Arithmetic and converting the 
           affine forms to intervals.
      
      "id" Butterfly enclosured obtained by computing {f}
           at the center of the sub-interval and the derivative
           {f'} over the whole sub-interval, both with IA,
      
      "aa" Parallelogram enclosures obtained by evaluating {f}
           with Affine Arithmetic. 
           
    The {suffix} may be "p" for a plot that shows the progress of the 
    zero-finding algorithm, and "s" for a plot that shows the final
    subdivision of {xd} into sub-intervals, and the sign of {f}
    in each sub-interval.
    
    The parameter {yd} is used only to define the vertical scale and
    range of the plot.
    
    The "true" graph of {f} is overlaid on all the plots.  The
    parameter {m} is the number of steps to use in this plot. */

#endif







