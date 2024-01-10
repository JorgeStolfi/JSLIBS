/* FOI test: plotting zeros of a 2-argument function (EPS version). */

#ifndef FOI2QE_H
#define FOI2QE_H

#include "foi.h"
#include "foifloat.h"
#include "interval.h"
#include <stdio.h>

void foi2qe_plots(
    char *filename1, char *filename2,
    char *title,
    Float f(Float x, Float y),
    Interval fv(Interval x, Interval y),
    FOIP ff(FOIP x, FOIP y),
    Interval xd,
    Interval yd,
    int n,
    int m
  );
  /* Generates two Encapsulated PostScript files with plots of the */
  /* curve f(x,y)=0, within the rectangle /xd/ by /yd/. */
  /* /fv/ and /ff/ should be Interval and FOI versions of /f/. */
  /* /n/ is the grid size used for testing /fv/ and /ff/. */
  /* /m/ is the grid size used for plotting the "true" zeros of /f/. */

#endif




