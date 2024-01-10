/* FOI test: plotting graph of a 1-argument function. */

#ifndef FGRAPH_H
#define FGRAPH_H

#include "foifloat.h"
#include "interval.h"
#include <stdio.h>

void fgraph_plot(
    FILE *psfile,
    Float f (Float x),
    Interval xd,
    Interval yd,
    int m
  );
  /* Draws the graph of $f(x)$ for $x \in xd$. */
  /* Assumes $*psfile$ is a Postscript output file,
     properly initialized by $plt0$ for a $xd \x yd$ plot. */
  /* The graph is approximated by $m$ straight segments. */

void fgraph_draw_axes(
    FILE *psfile,
    Interval xd,
    Interval yd
  );
  /* Draws the coordinates axes. */
  /* Assumes $*psfile$ is an output Postcript file,
     properly initialized by $plt0$ for a $xd \x yd$ plot. */
  
#endif
