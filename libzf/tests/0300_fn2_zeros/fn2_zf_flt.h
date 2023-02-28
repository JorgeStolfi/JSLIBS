/* Plotting the "true" zeros of a 2-argument function. */
/* Last edited on 2023-02-18 04:25:52 by stolfi */ 

#ifndef fn2_zf_flt_H
#define fn2_zf_flt_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <flt.h>
#include <ia.h>
#include <epswr.h>

void fn2_zf_flt_plot
  ( epswr_figure_t *eps,
    Float f (Float x, Float y),
    Interval xd,
    Interval yd,
    int32_t m
  );
  /* Plots "true" roots of {f(x,y) = 0}, in the rectangle {xd×yd}.
    Assumes {eps} is a properly initialized Encapsulated Postscript file. 
    The "true" zeros are found by evaluating {f} at the vertices and
    mid-cell points an {m×m} regular grid. */

#endif




