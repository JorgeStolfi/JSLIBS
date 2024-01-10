/* Plotting the "true" zeros of a 2-argument function. */
/* Last edited on 2007-12-26 21:29:49 by stolfi */ 

#ifndef fn2_zf_flt_H
#define fn2_zf_flt_H

#include <flt.h>
#include <ia.h>
#include <pswr.h>

#include <stdio.h>

void fn2_zf_flt_plot(
    PSStream *ps,
    Float f (Float x, Float y),
    Interval xd,
    Interval yd,
    int m
  );
  /* Plots "true" roots of $f(x,y) = 0$, in $xd \x yd$. */
  /* Assumes $psfile$ is a Postscript file, initialized by $plt0$. 
     The "true" zeros are found by evaluating $f$ at the vertices and
     mid-cell points an $m \x m$ regular grid. */

#endif




