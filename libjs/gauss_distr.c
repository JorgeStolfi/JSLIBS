/* See gauss_distr.h */
/* Last edited on 2019-12-05 20:03:22 by jstolfi */

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
    if (dev == INF) { return (x == +INF ? NAN : 0.0); }
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
    if (dev == INF) { return ((x0 == +INF) || (x1 == +INF) ? NAN : 0.0); }
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
    
    if ((z0 >= 0.0) && (z1 >= 0.0))
      { /* Use the difference of {erfc} because it is more accurate: */
        return (erfc(z0) - erfc(z1))/2;
      }
    else if ((z0 <= 0.0) && (z1 <= 0.0))
      { /* Use the difference of {erfc} reversed because it is more accurate: */
        return (erfc(-z1) - erfc(-z0))/2;
      }
    else
      { /* No point in using {erfc}: */
        /* !!! Should check for {t0} near {t1} and try to get more accurate !!! */
        return (erf(z1) - erf(z0))/2;
      }
  }

