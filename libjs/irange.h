#ifndef irange_H
#define irange_H

/* Intervals with {integer} endpoints. */
/* Last edited on 2007-05-09 08:39:30 by stolfi */ 
 
#include <stdint.h>

/* IRANGES

  An {irange_t} {I} is a pair of {int32_t}s, {LO(I) = I.end[BLO]} and
  {HI(I) = I.end{BHI}}. 
  
  When used to represent a range of real numbers, it usually means all
  real numbers {x} such that {LO(I) < x < HI(I)}. Depending on the
  context, the endpoints {LO(I)} and/or {HI(I)} may be included too.
  
  When used to represent a range of integers, it usually means {a..b}
  where {a} is eiher {I.end[0]} or {I.end[0]+1}, and {b} is either
  {I.end[1]} or {I.end[1]-1}, depending on the context. */
  
typedef struct irange_t { int32_t end[2]; } irange_t;
   /* An integer interval. */
     
#define ILO(x) ((x).end[0])
#define IHI(x) ((x).end[1])
  /* Handy macros for the two endpoints of an {irange_t}. */

typedef enum { irange_BLO = 0, irange_BHI = 1 } irange_side_t;
  /* Binary direction along an axis : {BLO} towards {-oo}, {BHI}
    towards {+oo}. Used e. g. to identify the two halves or endpoints
    of an irange. */

irange_t irange_join(irange_t *u, irange_t *v);
  /* Returns the smallest irange enclosing both {u} and {v},
  namely {(irange_t){min(LO(u),LO(v)),max(HI(u),HI(v))}}. */

irange_t irange_meet(irange_t *u, irange_t *v);
  /* Returns the intersection of {u} and {v},
  namely {(irange_t){max(LO(u),LO(v)),min(HI(u),HI(v))}}. */

void irange_widen(irange_t *r, int32_t margin);
  /* Widens {*r} by the specified {margin} on both sides.
    The {margin} may be negative. */

#endif
