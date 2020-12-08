/* See jsroot.h */
/* Last edited on 2008-10-05 18:31:21 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <jsroot.h>

double jsroot_bisect_secant(double xLo, double xHi, double (*f)(double x))
  {
    /* Values of {rel_terms(f)} at {xLo,xHi}: */
    double fLo = f(xLo);  if (fLo == 0) { return xLo; }
    double fHi = f(xHi);  if (fHi == 0) { return xHi; }
    demand(fLo*fHi <= 0, "the given interval does not bracket the root");
    double xdir = (xLo < xHi ? +1 : -1); /* Sign of mean derivative of {f} in interval. */
    double fdir = (fLo < fHi ? +1 : -1); /* Sign of mean derivative of {f} in interval. */
    while (TRUE)
      { assert(xdir*xLo < xdir*xHi);
        if (fabs(xHi - xLo) <= 1.0e-8) { break; }
        /* fprintf(stderr, "    x = [ %18.10f  %18.10f ]  f = [ %18.10f  %18.10f ]\n", xLo, xHi, fLo, fHi); */
        assert((fdir*fLo <= 0) && (fdir*fHi >= 0));

        /* Bissection step: */
        double xN = (xLo + xHi)/2;
        double fN = f(xN);  if (fN == 0) { return xN; }
        if (fdir*fN < 0)
          { xLo = xN; fLo = fN; }
        else
          { xHi = xN; fHi = fN; }

        assert(xdir*xLo < xdir*xHi);
        if (fabs(xHi - xLo) <= 1.0e-8) { break; }
        /* fprintf(stderr, "    f = [ %18.10f  %18.10f ]", xLo, xHi); */
        /* fprintf(stderr, "  a = [ %18.10f  %18.10f ]", fLo, fHi); */
        /* fprintf(stderr, "  g = [ %18.10f  %18.10f ]\n", fLo, fHi); */
        assert((fdir*fLo <= 0) && (fdir*fHi >= 0));

        /* Secant step: */
        double rG = (0 - fLo)/(fHi - fLo);
        double sG = 1 - rG;
        double xG = sG*xLo + rG*xHi;
        double fG = f(xG);  if (fG == 0) { return xG; }
        if (fdir*fG < 0)
          { xLo = xG; fLo = fG; }
        else
          { xHi = xG; fHi = fG; }
      }
    return (xHi + xLo)/2;
  }
