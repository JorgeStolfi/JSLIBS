#ifndef zf_H
#define zf_H

/* General interval-based zero finder for a 1-argument function. */
/* Last edited on 2023-03-18 10:42:34 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>
#include <flt.h>
#include <ia.h>
#include <ia_butfly.h>

typedef enum
  { zf_kind_undefined = 0,  /* F undefined in whole interval */
    zf_kind_positive = 1,   /* F positive in whole interval */
    zf_kind_negative = 2,   /* F negative in whole interval */
    zf_kind_root = 3,       /* F zero in whole interval, or interval is small */
    zf_kind_mixed = 4       /* Used internally */
  } zf_kind_t;

typedef void (zf_eval_func_t) (Interval *xr, Interval *yr, ia_butfly_t *a);
  /* Type of a procedure that can be passed as the {eval} argument
    of {zf_enum_zeros}, to evaluate the target function {F}. 
    
    Given an interval {*xr}, it should return an interval {*yr} that
    contains all values {F(x)} such that {x} is in {*xr} and {F(x)} is
    defined. The routine should also return a butterfly {*a} that also
    contains all the points {(x,F(x))}, for those values of {x}. Note
    that {*a} may be a trivial butterfly consisting of the box {xr × yr}.

    The routine {eval} may set {*yr} to an empty interval if the
    target function {F} is undefined in the whole interval {*xr}. */

typedef bool_t (zf_report_func_t)(Interval *xr, Interval *yr, zf_kind_t kind);
  /* Type of a procedure that can be passed as the {report} argument
    of {zf_enum_zeros}, to consume an interval {*xr} found by the 
    procedure. */

zf_kind_t zf_enum_zeros
  ( zf_eval_func_t *eval,
    zf_report_func_t *report,
    Interval xd,
    double epsilon,
    double delta
  );
  /* Finds the roots of a function {Y=F(X)} for {x} in {xd}, given 
    a procedures {eval} that bounds the graph of {F} on arbitrary
    sub-intervals of {xd}.
    
    More precisely, the {zf_enum_zeros} routine decomposes the interval
    {xr} into the union of zero or more intervals of non-zero length
    with pairwise disjoint interiors.  Each interval is classified as
    either "root", "positive", "negative", or "undefined".
    
    An interval {xr} is classified "positive" or "negative"
    if {F} is known to have that sign in the whole interval {xv}.  
    It is classified as "undefined" if the {eval} routine says so (see
    below).  Finally, it is classified as "root" if none
    of these conditions holds, and either {F} known to be zero in
    the whole {xv}, or the interval is small enough, that is
    
      { xv.hi - xv.lo  <=  max(delta, epsilon * max(|xv.lo|, |xv.hi|)) }

    The routine {zf_enum_zeros} enumerates the intervals of such a partition
    in increasing order.  This stream of intervals is passed though a
    filter that merges consecutive intervals of the same kind, and
    is then sent to the client via the {report} procedure.

    Each call to {report} is given the next (filtered) interval {xr},
    its kind {k}, and an interval {yr} that contains all the values
    {F(x)} where {x} is in {xr} and {F(x)} is defined. For "undefined"
    intevals, the range {yr} will be empty. For "positive" and
    "negative" intervals, the range {yr} may extend all the way to
    zero, even though the function {F} is guaranteed not to be zero in
    that interval.
    
    Note that since successive intervals share an endpoint, 
    between a "positive" and a "negative" interval there must be
    at least one "root" or "undefined" interval.

    If the {report} procedure returns TRUE, {zf_enum_zeros} exits
    immediately.

    In any case, {zf_enum_zeros} returns as a result the kind of the 
    interval that would follow the last reported interval; or
    {zf_kind_undefined} if the last reported interval 
    ends at {xd.hi}.
    
    The function {F} must be coded as a routine {eval} with the
    properties explained under {zf_eval_func_t}. The argument {xr}
    given to {eval} is always a sub-interval of {xd}.

    The {data} argument is a client-defined working storage and
    parameter area, that is passed to the {eval} and {report}
    routines.  */

/* HANDY ROUTINES */

zf_kind_t zf_classify_interval
  ( Interval *xr,
    Interval *yr,
    Float max_root_hi
  );
  /* Returns the kind of the interval {xr}, assuming that the
    values of {F} in it are contained in the interval {yr}.
    Note: if {yr} is empty, sets it to the standard empty interval
    [1 __ 0], and returns zf_kind_undefined. */
  
void zf_refine_interval(ia_trapez_t *tr, /*OUT:*/ Interval *xm, Interval *ym);
  /* Returns in {xm} a sub-interval of {xr} containing all roots of {F(X)} in
    {xr}, assuming that its graph is contained in the trapezoid whose
    altitude is {xr} and whose bases are {yxlo} at {xr.lo} and {yxhi} at
    {xr.hi}. also returns in {ym} an interval containing {F(X)}
    for all {X} in {xm}. */

#endif







