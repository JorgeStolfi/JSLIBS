/* aagraph.h -- Plots graph of unary function with affine arithmetic. */
/* Last edited on 2007-12-26 14:17:45 by stolfi */

#ifndef aagraph_H
#define aagraph_H

#include <flt.h>
#include <aa.h>
#include <ia.h>
#include <pswr.h>

void aagraph_plot_paralelograms
  ( PSStream *ps,
    AAP f (AAP x),
    Interval xd,
    Interval yd,
    int n
  );
  /* Writes to {ps} a PostScript plot of the graph of f(x),
    computed on {n} equal sub-intervals in {xd}.
    Specifically, converts each sub-interval {xi} into an affine form
    {xa}, with a new noise variable, and then plots the paralogram
    defined by the forms {(xa, f(xa))}.
    
    Assumes that {ps} has been properly initialized for a plot
    with the given coordinate ranges. */

void aagraph_plot_boxes
  ( PSStream *ps,
    AAP f (AAP x),
    Interval xd,
    Interval yd,
    int n
  );
  /* Writes to {ps} a PostScript plot of the graph of f(x),
    computed on {n} equal sub-intervals in {xd}. Specifically,
    converts each sub-interval {xi} into an affine form {xa}, with a
    new noise variable, and then plots the rectangle defined by the
    intervals {(xi, aa_range(f(xa)))}. Intervals are clipped to {yd}
    in the vertical direction.
    
    Assumes that {ps} has been properly initialized for a plot
    with the given coordinate ranges. */

void aagraph_fill_and_draw_2d_range 
  ( PSStream *ps, 
    AAP x, 
    AAP y, 
    double R, double G, double B
  );
  /* Writes to {ps} commands that draw the joint range of the
    affine forms {x} and {y}, filled with the given {R,G,B} color and
    outlined with the current pen. Assumes that {ps} has been
    properly initialized. */
 
#endif
