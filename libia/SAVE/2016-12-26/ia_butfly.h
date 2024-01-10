/* Tolls for "butterfly" enclosures of univariate functions. */
/* Last edited on 2007-12-26 15:05:44 by stolfi */
/* Created by Jorge Stolfi 93-01-13             */

#ifndef ia_butfly_H
#define ia_butfly_H

#include <flt.h>
#include <ia.h>
#include <ia_trapez.h>
#include <pswr.h>

typedef struct ia_butfly_t { ia_trapez_t tp[2]; } ia_butfly_t;
  /* A region of the plane consisting of two trapezoids
    {tp[0],tp[1]} with adjacent X ranges.  The reason for
    the name is the interval-slope arithmetic model where
    the region is usually narrower in the middle.
    
    The trapezoids should have {tp[0].x.hi == tp[1].x.lo}. Moreover,
    the intervals {tp[0].yxhi} and {tp[1].yxlo} should overlap. */

ia_butfly_t ia_butfly_from_box(Interval *xr, Interval *yr);
  /* Converts an IA graph enclosure (the rectangle {xr × yr})
    into a {ia_butfly_t}. The low trapezoid of the butterfly 
    is a vertical segment at {xr.lo}, and the high trapezoid is the
    rectangle {xr × yr}. */

ia_butfly_t ia_butfly_from_ia_diff(Interval *xr, Float xm, Interval *yxmr, Interval *dyr);
  /* Converts a value-slope graph enclosure into a {ia_butfly_t}. 
  
    Let {F} be any continuous and differentiable function {F} such
    that {F(xm)} is contained in {yxmr} and the derivative {F'(x)} is
    contained in {dyr}, for any {x} in {xr}. The returned butterfly
    {a} contains all the points {(x,F(x))} such that {x} is in {xr}.
    
    This procedure only works if {F(x)} is defined over all 
    points of {xr}, or if {yxmr} is {ia_full()}, or if {dyr} is 
    {ia_full()}. */

ia_butfly_t ia_butfly_from_trapez(ia_trapez_t *tp);
  /* Converts a trapezoidal graph enclosure {*tp} into a
    {ia_butfly_t}. The resulting butterfly has an empty low trapezoid,
    and its high trapezoid is guaranteed to contain the part of {*tp}
    whose X-projetion is {xr}. */

void ia_butfly_print (FILE *wr, ia_butfly_t *bt, char *sep);
  /* Prints {bt} to {wr}.  Each trapezoid is printed
    with {ia_trapez_print}, separated by the string {sep}. */
    
void ia_butfly_draw(PSStream *ps, Interval *yr, ia_butfly_t *bt);
  /* Draws the butterfly {bt} to Postscript file {ps}, as
    a pair of trapezoids. If the Y range of either trapezoid
    is full, plots a box with Y range {yr} instead. */

#endif
