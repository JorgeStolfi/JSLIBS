/* Primitives for the Boareto 113 side porch drawings */
/* Last edited on 2024-12-05 10:17:06 by stolfi */

#ifndef boap_elems_H
#define boap_elems_H

#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <sign.h>
#include <epswr.h>

/* All dimensions are in mm.

  The World coordinate axes are: {X} horizontal towards the observer,
  {Y} horizontal, left to right seen from {+X}, {Z} up. The Plot
  coordinate axes are: {x} horizontal, left to right, {y} vertical,
  up. */
  
epswr_figure_t *boap_new_figure
  ( char *pname, 
    char *subname, 
    double xMin, double xMax, 
    double yMin, double yMax, 
    bool_t landscape
  );
  /* Starts a new EPS figure file named "out/porch_{pname}_{subname}.eps". The figure bounding box
    will be approximately that of a Letter page, upright or landscape as
    requested. The plot area will be the bounding box minus a small
    margin, with aspect ratio 10:13 or 13:10. The client to device
    mapping will be set to that the rectangle {[xMin _ xMax] × [yMin _ yMax]}
    will fit in the plot window. */
    
#endif
