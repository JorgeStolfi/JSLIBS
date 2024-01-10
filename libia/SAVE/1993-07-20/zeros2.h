/* Plotting the "true" zeros of a 2-argument function. */

#ifndef ZEROS2_H
#define ZEROS2_H

#include "foifloat.h"
#include "interval.h"
#include <stdio.h>

void zeros2_plot(
    FILE *psfile,
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




