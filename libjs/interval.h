#ifndef interval_H
#define interval_H

/* Intervals with {double} endpoints. */
/* Last edited on 2012-01-04 00:23:39 by stolfi */ 

/* Should be merged with (or replaced by) {ia.h}. */
 
/* INTERVALS
  
  An {interval_t} {I} is a pair of {double}s, the /inferior/ or /low
  endpoint/ {LO(I) = I.end[0]} and the /superior/ or /high endpoint/
  {HI(I) = I.end[1]}.  The endpoints may be {±INF} but not {NAN}.
  
  The procedures in this interface generally assume that the intervals
  are closed, that is, that they represent the set of all extended
  real numbers {z} (all the real numbers of mathematics, plus the
  infinities {-INF} and {+INF}) such that {LO(X) <= z <= HI(X)}.
  
  Clients that assume other conventions for the endpoints (open or
  half-open) or that interpret intervals as sets of other numbers
  (finite reals, rationals, floats, etc.) should check carefully the
  procedure specs.  Thus, for example, the function {interval_is_empty}
  below returns TRUE iff the /closed/ interval is empty, that is,
  {LO(X) > HI(X)}; whereas an open or half-open interval is empty
  iff {LO(X) >= HI(X)}. */
  
#include <jsmath.h>

typedef struct interval_t { double end[2]; } interval_t;
   /* A real interval. */
     
#define LO(X) ((X).end[0])
#define HI(X) ((X).end[1])
  /* Handy macros for the two endpoints of an interval. */

typedef enum { interval_BLO = 0, interval_BHI = 1 } interval_side_t;
  /* Binary direction along an axis : {BLO} towards {-oo}, {BHI}
    towards {+oo}. Used e. g. to identify the two halves or endpoints
    of an interval. */

double interval_is_empty(interval_t *X);
  /* Returns TRUE iff the (closed) interval {X} is empty, that is,
    iff {LO(X) > HI(X)}. */

double interval_is_full(interval_t *X);
  /* Returns TRUE iff the (closed) interval {X} contains all extended reals, that is,
    iff {LO(X) == -INF} and {HI(X) == +INF}. */

double interval_is_finite(interval_t *X);
  /* Returns TRUE iff the (closed) interval {X} contains only finite reals
    (or is empty), that is, iff {LO(X) > -INF} and {HI(X) < +INF}. */

double interval_is_trivial(interval_t *X);
  /* Returns TRUE iff the (closed) interval {X} contains only one extended
    real number, that is, iff {LO(X) == HI(X)}. */

void interval_mid_rad (interval_t *X, double *mid, double *rad);
  /* If the closed interval {X} is not empty, returns the approximate
    midlpoint {mid} of {X} (which is guaranteed to be in {X}), and a
    non-negative value {rad} such that {X} is tightly contained in
    {[mid-rad _ mid+rad]}, that is, {mid-rad <= LO(X)} and {mid+rad >=
    HI(X)}. If the closed interval {X} is empty, returns {mid = 0, rad
    = -INF}.
    
    In particular, if {X} is full, returns {mid = 0, rad = +INF}. If
    {X} is trivial returns {mid = LO(X) = HI(X), rad = 0}. If {X} is
    not full but infinite, returns {mid = -INF} or {mid = +INF}, and
    {rad = +INF}. Finally, if {X} is finite, non-empty, and
    non-trivial, returns a finite {mid} inside {X}, and a strictly
    positive (perhaps infinite) {rad}. */

double interval_mid (interval_t *X);
  /* If the closed interval {X} is finite and not empty, 
    returns the approximate midpoint of {X}, guaranteed to be 
    finite and inside {X}.  If {X} is empty or full, returns 0.
    Otherwise, if {X} is finite.   */

double interval_rad (interval_t *X);
  /* If {X} is not empty, the radius of {X} from its midpoint {m =
    interval_mid(X)}. I.e. a value {r} such that {[m-r _ m+r]}
    contains {X}. Negative if {X} is empty. Finite as long as {X} is
    finite. */

interval_t interval_from_mid_rad (double mid, double rad);
  /* If {rad} is negative, returns an empty interval. If {rad} is
    {+INF}, returns the full interval {-INF,+INF} (even if {mid} is
    {±INF}. Otherwise returns an interval that includes {[mid-rad,
    mid+rad]}. If either endpoint overflows, sets it to infinity.*/

double interval_width (interval_t *X);
  /* The width of {X}, i.e. {HI(X) - LO(X)}, rounded up. Returns
    {+INF} if the subtraction overflows. */

interval_t interval_split(interval_t *X, interval_side_t dir);
  /* Returns the lower or upper half of the closed interval {X},
    depending on {dir}. The splitting point is {interval_mid(X)}, and
    is included in both halves. */

interval_t interval_join(interval_t *X, interval_t *Y);
  /* Returns the smallest interval enclosing both {X} and {Y}. */

interval_t interval_meet(interval_t *X, interval_t *Y);
  /* Returns the intersection of {X} and {Y}. */

void interval_widen(interval_t *X, double margin);
  /* Widens {*X} by the specified {margin} on both sides.
    The {margin} may be negative, in which case the interval
    is narrowed (and may become empty). */

void interval_adjust_ratio(interval_t *X, interval_t *Y, double tx, double ty);
  /* Widens either {*X} or {*Y}, as needed, to ensure that 
    their widths are in the ratio {tx:ty} (except from roundoff
    errors). */

double interval_project(interval_t *X, double y);
  /* Returns the value {z} in the closed interval {*X} which is 
    closest to {y} (which is {y} itself if {y} is inside {*X}, 
    otherwise it is the endpoint of {X} nearest to {y}). 
    Fails if {*X} is empty. */

#endif
