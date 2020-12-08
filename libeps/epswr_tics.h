/* Plots tic marks and labels along specified axes. */
/* Last edited on 2009-08-24 22:31:46 by stolfi */

#ifndef epswr_tics_H
#define epswr_tics_H

#include <epswr.h>
  
void epswr_tics
  ( epswr_figure_t *epsf, 
    double ox,
    double oy,
    double dx,
    double dy,
    double lo, 
    double hi,
    int n,
    char *fmt,
    double ticSize,
    double align 
  );
  /* Draws {n+1} coordinate tics along the line that goes
    through the point {o=(ox,oy)} and is parallel to the vector {d=(dx,dy)}.
    The first tic will be at distance {lo} from {o}, and the last one will
    be at distance {hi}; both may be negative.  If {fmt} is not null, also writes the
    corresponding coordinate values, with that format. The tics have
    length {ticSize} (in mm) and extend from {-align*ticSize} to
    {(1-align)*ticSize} in the direction that is 90 degrees counterclockwise from {d}.
    {ticsize} may be negative in which case. */

#endif
