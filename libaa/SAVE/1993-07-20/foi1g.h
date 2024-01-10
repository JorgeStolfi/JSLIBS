/* FOI test: plotting graph of a 1-argument function. */

#ifndef FOI1G_H
#define FOI1G_H

#include "foi.h"
#include "foifloat.h"
#include "interval.h"
#include <stdio.h>

void foi1g_plots(
    char *filename,
    char *title,
    Float f(Float x),
    Interval fv(Interval x),
    FOIP ff(FOIP x),
    Interval xd,
    Interval yd,
    int n,
    int m
  );
  /* Generates a PostScript file with plots of the graph of f(x). */
  /* /fv/ and /ff/ should be Interval and FOI versions of /f/. */
  /* /xd/ and /yd/ are the plot coordinate ranges. */
  /* /n/ is the number of intervals into which /xd/ is to be divided. */
  /* /m/ is the number of steps for plotting the "true" graph of /f/. */

#endif
