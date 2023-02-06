/* Tolls for trapezoidal enclosures of univariate functions. */
/* Last edited on 2023-02-03 22:16:38 by stolfi */
/* Created by Jorge Stolfi 93-01-13             */

#ifndef ia_trapez_H
#define ia_trapez_H

#include <flt.h>
#include <ia.h>

typedef struct ia_trapez_t { Interval x, yxlo, yxhi; } ia_trapez_t;
  /* A trapezoidal region of the plane {R^2}.  The projection 
    on the X-axis is the interval {x}.  The bases are vertical 
    segments, with abscissas {x.lo} (spanning the Y range {yxlo})
    and {x.hi} (spanning the Y range {yxhi}. */

ia_trapez_t ia_trapez_from_box(Interval *xr, Interval *yr);
  /* Repackages an IA graph enclosure (the rectangle {xr × yr})
    into a {ia_trapez_t}. */

ia_trapez_t ia_trapez_from_ia_diff(Interval *xr, Interval *yr, Interval *dyr, int which);
  /* Converts a value-slope graph enclosure into a {ia_trapez_t}. 
    
    The trapezoid has X projection equal to {xr} and two vertical sides.
    One of these sides is the Y interval {yr}: the ``low'' side
    (at {x.lo}) if {which} is 0, of the ``high'' side if {which} is 1.
    The opposide side is computed so as to enclose all straight lines
    that cut the first side and have slope contained in {dyr}. 

    Let {F} be any continuous and differentiable function {F} such
    that {F(x.lo)} or {F(x.hi)} is contained in {yr} and the
    derivative {F'(x)} is contained in {dyr}, for any {x} in {xr}. The
    returned trapezoid {tr} contains all the points {(x,F(x))} such
    that {x} is in {xr}.
    
    This procedure only works if {F(x)} is defined over all points of
    {xr}, or if {yr} is {ia_full()}, or if {dyr} is {ia_full()}. */

ia_trapez_t ia_trapez_clip(Interval *xr, ia_trapez_t *tp);
  /* Computes a trapezoidal enclosure for the intersection of
    trapezoid {tp} and the vertical band whose X-projetion is {xr}. */

void ia_trapez_print(FILE *wr, ia_trapez_t *tr);
  /* Prints trapezoid {tr} to file {wr}, in the format
    "[{x.lo}: {yxlo} __ {xb}: {yxhi}]", where {yxlo} and {yxhi}
    are printed with {ia_print}. */
    
void ia_trapez_fill(PSStream *ps, Interval *yr, ia_trapez_t *tr, float cr, float cg, float cb);
void ia_trapez_draw(PSStream *ps, Interval *yr, ia_trapez_t *tr);
  /* These procedures plot a trapezoid {tr} into a Postscript file
    {ps}. The procedure {ia_trapez_fill} fills the trapezoid with
    color {cr,cg,cb} in the RGB model, with each component in [0_1].
    The procedure {ia_trapez_draw} draws the trapezoid's outline with
    the current pen. If the Y range of the trapezoid is full, plots
    the box {tr->x × yr} instead. */

#endif
