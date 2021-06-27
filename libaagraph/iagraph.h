/* Graphing 1-argument functions with interval arithmetic. */
/* Last edited on 2021-06-26 02:28:10 by jstolfi */

#ifndef iagraph_H
#define iagraph_H

#include <flt.h>
#include <ia.h>
#include <ia_butfly.h>
#include <epswr.h>

void iagraph_plot_boxes
  ( epswr_figure_t *fig,
    Interval f (Interval x),
    Interval xd,
    Interval yd,
    int n
  );
  /* Writes to {ps} a box enclosure of the graph of a function
    {F(x)}, computed on {n} equal sub-intervals in {xd}. Specifically, for
    each sub-interval {xi}, plots the box {(xi,f(xi))}. The Y intervals
    are clipped to {yd}.
    
    Assumes that {ps} has been properly initialized for a plot
    with the given coordinate ranges. */

void iagraph_plot_butterflies
  ( epswr_figure_t *fig,
    Interval f (Interval x),
    Interval df (Interval x),
    Interval xd,
    Interval yd,
    int n
  );
  /* Writes to {ps} an interval-slope enclosure of the graph of a
    function {F(x)}, computed on {n} equal sub-intervals in {xd}.
    Specifically, for each sub-interval {xi}, plots the
    butterfly-shaped region implied by the value of {F} on the
    midpoint {xc} of the interval {xi} (as given by {f(xc)}), and by
    the range of the derivative of {F} over {xi} (as given by
    {df(xi)}). The butterflies are clipped in the Y drection to {yd}.
    
    Assumes that {ps} has been properly initialized for a plot
    with the given coordinate ranges. */

/* HANDY TOOLS */

void iagraph_compute_butterfly
  ( Interval xv, 
    Float xm, 
    Interval f (Interval x),
    Interval df (Interval x),
    ia_butfly_t *bt
  );
  /* Computes a butterfly-shaped region of the XY plane
    that encloses the graph of {F(x)} for all {x} in the 
    interval {xv}. 
    
    The region consists of two trapezoids with vertical bases,
    spanning interval {xv} along the X direction. The first trapezoid spans the
    X interval {[xv.lo _ xm]}, the second one spans {[xm _ xv.hi]}.
      
    The middle abscissa {xm} must belong to {xv}. If {xm==xv.lo} or
    {xm==xv.hi}, one of the trapezoids degenerates to a vertical
    segment coincident with the base of the other trapezoid. */

/* PLOTTING ENCLOSURES

  These procedures paint the specified enclosure with 
  color {R,G,B} and outline it with the current pen.
  
  They should clip the given shape against the `domain' box 
  {xd × yd}, but that feature is not fully implemented yet.
  !!! FIX IT !!! */

void iagraph_fill_and_draw_box
  ( epswr_figure_t *fig, 
    Interval xv, 
    Interval yv, 
    Interval xd,
    Interval yd,
    double R, double G, double B
  );
  /* Plots to {ps} the rectangular box {xv × yv}. */

void iagraph_fill_and_draw_trapezoid
  ( epswr_figure_t *fig, 
    ia_trapez_t *tp, 
    Interval xd,
    Interval yd,
    double R, double G, double B
  );
  /* Plots to {ps} a trapezoid {tp}. */
    
void iagraph_fill_and_draw_butterfly
  ( epswr_figure_t *fig, 
    ia_butfly_t *bt, 
    Interval xd,
    Interval yd,
    double R, double G, double B
  );
  /* Plots to {ps} a butterfly-shaped region {bt}.
    See {iagraph_compute_butterfly} above. */

#endif
