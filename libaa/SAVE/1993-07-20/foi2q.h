/* FOI test: plotting zeros of a 2-argument function. */

#ifndef FOI2Q_H
#define FOI2Q_H

#include "foi.h"
#include "foifloat.h"
#include "interval.h"
#include <stdio.h>

void foi2q_plots(
    char *filename,
    char *title,
    Float f(Float x, Float y),
    Interval fv(Interval x, Interval y),
    FOIP ff(FOIP x, FOIP y),
    Interval xd,
    Interval yd,
    int n,
    int m
  );
  /* Generates a PostScript file with plots of the curve f(x,y)=0, */
  /* within the rectangle /xd/ by /yd/. */
  /* /fv/ and /ff/ should be Interval and FOI versions of /f/. */
  /* /n/ is the grid size used for testing /fv/ and /ff/. */
  /* /m/ is the grid size used for plotting the "true" zeros of /f/. */

#endif
