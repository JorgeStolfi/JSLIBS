/* Plotting the graph of a 1-argument function. */
/* Last edited on 2021-06-26 02:27:26 by jstolfi */

#ifndef fltgraph_H
#define fltgraph_H

#include <flt.h>
#include <ia.h>
#include <epswr.h>

/* These procedures assume that {*fig} is a {epswr_figure_t} object, as defined
  in {epswr.h}, properly initialized for a plot with the given coordinate 
  ranges. */

void fltgraph_plot
  ( epswr_figure_t *fig,
    Float f (Float x),
    Interval xd,
    Interval yd,
    int m
  );
  /* Draws the graph of {f(x)} for {x \in xd},
    approximated by {m} straight segments. */

void fltgraph_draw_axes
  ( epswr_figure_t *fig,
    Interval xd,
    Interval yd
  );
  /* Draws the parts of the coordinates axes that lie
    in the specified rectangle. */
  
void fltgraph_draw_tics
  ( epswr_figure_t *fig,
    epswr_axis_t axis,
    Float lo, Float hi,
    int n,
    double ticsz,
    char *labfmt,
    double labalign,
    Interval clip
  );
  /* Draws {n+1} coordinate tics on the given axis, from {lo} to
    {hi}. The tics measure {ticsz} millimeters irrespective of the
    plot scale. If {labfmt} is not null, also writes the coordinate
    value, with that format. The parameter {labalign} is the relative
    point of the label that should be placed on the axis; i.e.,
    {labalign < 0} places the label on the positive side of the axis,
    {labalign > 1.0 places it on the negative side, and {labalign =
    0.5} places the label centered on the axis. Supresses tics that
    are outside the {clip} interval. */

#endif
