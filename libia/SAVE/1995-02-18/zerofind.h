/* General zero finder for a 1-argument function. */

#ifndef ZEROFIND_H
#define ZEROFIND_H

#include <ia.h>
#include <flt.h>
#include <stdio.h>

typedef struct {Interval a, b; } IntervalPair;

typedef enum {
    zf_type_undefined = 0,  /* F undefined in whole interval */
    zf_type_positive = 1,   /* F positive in whole interval */
    zf_type_negative = 2,   /* F negative in whole interval */
    zf_type_root = 3,       /* F zero in whole interval, or interval is small */
    zf_type_complex = 4     /* Used internally */
  } zf_type;

zf_type zerofind(
    IntervalPair f (Interval *xv, void *data),
    int report (Interval *xv, Interval *yv, zf_type tv, void *data),
    void *data,
    Interval xd,
    Float epsilon,
    Float delta
  );
  /* 
    Finds the roots of a function $y=F(x)$ for $x$ in $xd$, given a
    procedure $f$ that bounds the graph of $F$ on arbitrary
    sub-intervals of $xd$.
    
    More precisely, the $zerofind$ routine decomposes the interval
    $xd$ into the union of zero or more intervals of non-zero length
    with pairwise disjoint interiors.  Each interval is classified as
    either "root", "positive", "negative", or "undefined".
    
    An interval $xv$ is classified "positive" or "negative"
    if $F$ is known to have that sign in the whole interval $xv$.  
    It is classified as "undefined" if the $f$ routine says so (see
    below).  Finally, it is classified as "root" if either
    $F$ known to be zero in the whole $xv$, or 
    the interval is small enough, that is
    
      xv.hi - xv.lo  <=  max(delta, epsilon * max(|xv.lo|, |xv.hi|))

    The $zerofind$ routine computes the intervals of such a partition
    in increasing order.  This stream of intervals is passed though a
    filter that merges consecutive intervals of the same type, and
    is then sent to the client via the $report$ procedure.

    Each call to $report$ is given the next (filtered) interval $xv$,
    its type $tv$, and a range $yv$ containing the corresponding $F$
    values.  For "undefined" intevals, the range $yv$ will be
    empty. For "positive" and "negative" intervals, the range $yv$ may
    extend all the way to zero, even though the function $F$ is
    guaranteed not to be zero in that interval.
    
    Note that since successive intervals share an endpoint, 
    between a "positive" and a "negative" interval there must be
    at least one "root" or "undefined" interval.

    If the $report$ procedure returns TRUE, $zerofind$ exits
    immediately.

    In any case, $zerofind$ returns as a result the type of the 
    interval that would follow the last reported interval; or
    zf_type_undefined if the last reported interval 
    ends at $xd.hi$.
    
    The function $F$ must be coded as a routine $f$ that, given an
    interval $xv$ for $x$, computes two intervals $yp.a$ and $yp.b$ 
    such that the graph of $F$ in $xv$ is contained in the 
    trapezoid with vertices 

      (xv.lo, yp.a.hi) (xv.lo, yp.a.lo), 
      (xv.hi, yp.b.lo), (xv.hi, yp.b.hi) 

    The routine $f$ may return a pair of empty intervals if the
    function $F$ is to be considered undefined in the interval $xv$.

    The $data$ argument is a client-defined working storage and
    parameter area, that is passed to the $f$ and $report$ routines.
    */

#endif







