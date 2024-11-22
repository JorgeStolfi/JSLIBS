/* See gauss_distr.h */
/* Last edited on 2024-11-15 19:12:37 by stolfi */

#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <values.h>
#include <assert.h>

#include <bool.h>
#include <jsmath.h>
#include <affirm.h>
#include <gauss_bell.h>
#include <gauss_distr.h>

double gauss_distr_PDF(double x, double avg, double dev)
  {
    dev = fabs(dev);
    if (dev == INF) { return 0.0; }
    if (dev == 0.0) { return (x == avg ? (fabs(avg) == INF ? NAN : +INF) : 0.0); }
    double z = fabs(x - avg)/dev;
    if (z >= gauss_distr_HUGE_ARG)
      { return 0.0; }
    else 
      { double A = dev*sqrt(2*M_PI);
        if (z < gauss_distr_TINY_ARG)
          { return 1.0/A; }
        else
          { return exp(-z*z/2)/A; }
      }
  }
  
double gauss_distr_CDF(double x, double avg, double dev)
  { dev = fabs(dev);
    if (dev == INF) { return (x == +INF ? 1.0 : 0.0); }
    if (dev == 0.0) { return (x < avg ? 0.0 : (x > avg ? 1.0 : NAN)); }
    if (x == avg)
      { return (fabs(avg) == INF ? NAN : 0.5); }
    else 
      { /* There seems to be no point in checking for too big args, {erf} does that itself: */
        double dev2 = dev * M_SQRT2;
        double z = (x - avg)/dev2; 
        return erfc(-z)/2;  /* More accurate for {x < avg} than {(1 + erf(t))/2}. */
      }
  }

double gauss_distr_integral(double x0, double x1, double avg, double dev)
  { /* There seems to be no point in checking for too big args, {erf} does that itself: */
    dev = fabs(dev);
    if (dev == INF) 
      { /* PDF is flat at 0, but integral is 1: */
        if (x0 == x1)
          { return 0.0; }
        else if ((x0 == -INF) && (x1 == +INF))
          { return +1.0; }
        else if ((x0 == +INF) && (x1 == -INF))
          { return -1.0; }
        else if ((fabs(x0) == INF) || (fabs(x1) == INF))
          { return NAN; }
        else
          { return 0.0; }
      }
    if (dev == 0.0) 
      { /* Dirac pulse: */
        if ((x0 < avg) && (x1 > avg))
          { return +1.0; }
        else if ((x0 > avg) && (x1 < avg))
          { return -1.0; }
        else if ((x0 == avg) || (x1 == avg))
          { return NAN; }
        else
          { return 0.0; }
      }
    double dev2 = dev * M_SQRT2;
    double z0 = (x0 - avg)/dev2;
    double z1 = (x1 - avg)/dev2;
    
    /* If {z>=1}, {1-erfc(z)} is more accurate than {erf(z)}. */
    double a0, c0; /* Such that {erf(z0)} is {a0+c0}, where {a0} is 0 or {±1}. */
    if (fabs(z0) <= 1)
      { a0 = 00; c0 = erf(z0); }
    else if (z0 > 0)
      { a0 = +1; c0 = -erfc(z0); }
    else
      { a0 = -1; c0 = +erfc(-z0); }

    double a1, c1; /* Such that {erf(z1)} is {a1+c1}, where {a1} is 0 or {±1}. */
    if (fabs(z1) <= 1)
      { a1 = 00; c1 = erf(z1); }
    else if (z1 > 0)
      { a1 = +1; c1 = -erfc(z1); }
    else
      { a1 = -1; c1 = +erfc(-z1); }
      
    return ((a1 - a0) + (c1 - c0))/2;
  }

