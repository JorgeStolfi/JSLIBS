/* See gauss_bell.h */
/* Last edited on 2019-12-05 10:27:06 by jstolfi */

#define _GNU_SOURCE
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <values.h>
#include <assert.h>

#include <bool.h>
#include <jsmath.h>
#include <affirm.h>
#include <gauss_bell.h>

double gauss_bell_eval(double x, double avg, double dev)
  {
    double y = fabs(x - avg);
    dev = fabs(dev);
    if (dev == INF)
      { return 1.0; }
    else if (dev == 0.0)
      { return (y == 0 ? 1.0 : 0.0); }
    else 
      { double z = y/dev;
        if (z >= gauss_bell_HUGE_ARG)
          { return 0.0; }
        else if (z < gauss_bell_TINY_ARG)
          { return 1.0; }
        else
          { return exp(-z*z/2); }
      }
  }
