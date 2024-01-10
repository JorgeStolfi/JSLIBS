/* General zero finder for a 1-argument function. */

#ifndef ZEROFIND_H
#define ZEROFIND_H

#include "ia.h"
#include <flt.h>
#include <stdio.h>

typedef struct {Interval a, b; } IntervalPair;

void zerofind(
    IntervalPair f (Interval xv, void *closure),
    void report (Interval xv, IntervalPair yv, void *closure),
    void *closure,
    Interval xd,
    Float epsilon,
    Float delta
  );
  /* 
    Finds the roots of a function $y=F(x)$ for $x$ in $xd$.
    
    The function must be coded as a routine $f$ that, given an
    interval $xv$ for $x$, computes two intervals $ya$ and $yb$ 
    such that the graph of $F$ in $xv$ is contained in the 
    trapezoid with vertices 

      (xv.lo, ya.hi) (xv.lo, ya.lo), 
      (xv.hi, yb.lo), (xv.hi, yb.hi) 

    An interval $xv$ is considered a root if $ya$ and $yb$ are
    both [0,0], or if
    
      xv.hi - xv.lo  <=  max(delta, epsilon * max(|xv.lo|, |xv.hi|))

    Consecutive root intervals are merged.  The $report$ procedure is
    called to process the resulting disjoint root intervals, in order of
    increasing $x$.

    The routine $f$ may return a pair of empty intervals if the
    function $F$ is undefined in the interval $xv$, or if that
    interval is uninteresting and should to be dropped from further
    consideration for any rason.  

    The $closure$ argument is a client-defined working storage and
    parameter area for the $f$ and $report$ routines.
  */

#endif







