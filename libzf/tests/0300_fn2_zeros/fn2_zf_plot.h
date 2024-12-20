/* Last edited on 2024-12-05 10:41:25 by stolfi */ 

#ifndef fn2_zf_plot_H
#define fn2_zf_plot_H

#include <stdio.h>
#include <stdint.h>

#include <ia.h>
#include <flt.h>
#include <fn2_functions.h>
#include <epswr.h>

epswr_figure_t *fn2_zf_plot_new_figure
  ( char *prefix, 
    char *arith_tag, 
    Interval xd,
    Interval yd,
    int32_t ncap,
    char *title, 
    char *parm_string, 
    char *arith_title
  );
  /* Opens a new EPS file called "out/{prefix}_{arith_tag}.eps"
    dimensioned to plot the rectangle {xd×yd}.
    Writes {title}, {parm_string}, and {arith_title} in the caption. */
    
void fn2_zf_plot_stat_caption
  ( epswr_figure_t *eps, 
    char *fmt, 
    int32_t count,
    char *arith
  );
  /* Writes a caption line containing the integer {count} and the string {arith}
    textified with the given format string {fmt}. Also writes that line to {stderr}.
    
    if {arith} is not {NULL}, the format {fmt} must contain a "%d" spec
    and a "%s" spec, in that order. If {arith} is {NULL}, that item is not
    printed, and in that case the {fmt} must have a single "%d" spec. */ 

void fn2_zf_plot_end_figure
  ( epswr_figure_t *eps,
    eval_fp_t *eval_fp,
    Interval xd,
    Interval yd,
    int32_t m, 
    int32_t cols, 
    int32_t rows
  );
  /* Finishes the figure {eps}. Plots the true zeros using {eval_fp},
    draws the grid lines for an array of {rows} by {cols} cells,
    draws the plot frame, and closes the EPS file. */
   
char *fn2_zf_plot_format_parms(
    Interval xd,
    Interval yd,
    int32_t n
  );
  /* Formats the arguments into a string */
 
#endif
