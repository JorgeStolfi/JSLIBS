#define PROG_NAME "test_gauss_distr"
#define PROG_DESC "test of {gauss_distr.h}"
#define PROG_VERS "1.0"

/* Last edited on 2019-12-05 20:24:31 by jstolfi */ 
/* Created on 2012-03-04 by J. Stolfi, UNICAMP */

#define test_hermite3_COPYRIGHT \
  "Copyright © 2017  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <gauss_distr.h>

#include <bool.h>
#include <jsmath.h>
#include <jsfile.h>
#include <affirm.h>

void test_basics(void);
  /* Runs tests on the functions of the interface {gauss_distr.h}. */

void test_all(double z);
  /* Tests the PDF and CDF functions for argument 
    {x = avg + dev*z}, and the two-point integral functions for arguments 
    {x} and {x1 = avg + dev*z1}, for various values of {z1}, {avg}, and {dev}. */

void test_func(int32_t which, double z, double z1, double v_pre);
  /* Tests a Gaussian distr function at arguments 
    {x = avg + dev*z} and {x1 = avg + dev*z1}, for various values of 
    {avg} and {dev}. Compares the result with {v_pre}, in log scale.
  
    If {which} is 0, tests {gauss_distr_PDF(x, avg, dev)}.
    
    If {which} is 1, tests {gauss_distr_CDF(x, avg, dev)},
    
    If {which} is 2, tests {gauss_distr_integral(x, x1, avg, dev)}.
    
    The {dev} is chosen with various signs and may be infinite or zero.
    When testing {dev=0}, picks "random" {x} and {x1} away from {avg}. 
    In that case, {v_pre}  is ignored and recomputed internally.  */

void test_check
  ( char *fname, 
    double x, 
    double x1, 
    double avg, 
    double dev, 
    double v_cmp,
    double v_pre
  );
  /* Checks the computed result {v_cmp} of some function
    of argument {x} (or {x} and {x1}, if {x1} is not {NAN})
    with the expected value {v_pre}. The string {fname} 
    identifies which function was called. */

double test_arg(int32_t k, int32_t N, double vmin, double vmax);
  /* Returns a test value from a series. If {k} is in {-N .. +N},
    returns {(-1)^k*vmin*(vmax/vmin)^k/N}. If {k} is in {-N-3}, {-N-2}, or {-N-1},
    returns {-INF}, {0}, or {+INF}, respectively. */

int main(int argn, char **argv);

/* IMPLEMENTATIONS */

int main (int argc, char **argv)
  {
    test_basics();
    return 0;
  }

double test_arg(int32_t k, int32_t N, double vmin, double vmax)
  { double v;
    if (k < 0)
      { v = (2.0 - (k + N)) * INF; }
    else
      { double s = 1.0 - 2*(k % 2);
        double fr = ((double)k)/((double)N);
        v = s * vmin * pow(vmax/vmin, fr);
      }
    return v;
  }

void test_basics(void)
  {
    int32_t Nz = 10;
    double zmin = 1.0e-6, zmax = 1.0e+6;
    for (int32_t kz = -3; kz <= Nz; kz++)
      { double z = test_arg(kz, Nz, zmin, zmax);
        test_all(z);
      }
   }

void test_all(double z)
  { double v_pre_PDF = exp(-z*z/2)/sqrt(2.0*M_PI);
    test_func(0, z, NAN, v_pre_PDF);
    double v_pre_CDF = 0.5*(1+ erf(z/M_SQRT2));
    test_func(1, z, NAN, v_pre_CDF);
    int32_t Ndz1 = 5;
    double dz1min = 1.0e-6, dz1max = 1.0e+6;
    for (int32_t kdz1 = -3; kdz1 <= Ndz1; kdz1++)
      { double dz1 = test_arg(kdz1, Ndz1, dz1min, dz1max);
        double z1 = z + dz1;
        double v_pre_integral = 0.5*erf(z1/M_SQRT2) - 0.5*erf(z/M_SQRT2);
        test_func(2, z, z1, v_pre_integral);
      }
  }
   
void test_func(int32_t which, double z, double z1, double v_pre)
  { /* Range of deviations to consider: */
    double dmin = 1.0e-6, dmax = 1.0e+6;
    int Nd = 3;
    double avg = 20.0;
    char *fname = NULL;
    for (int32_t kd = - Nd - 3; kd <= + Nd; kd++)
      { /* Choose a standard deviation {dev}: */
        double dev = test_arg(kd, Nd, dmin, dmax);
        double x, x1; /* Argument(s) of generic function. */
        double v_pre1 = NAN; /* Predicted value for this {dev}. */
        if (dev != 0.0)
          { /* PDF is bell-shaped: */
            /* Pick {x,x1} that correspond to the standard {z,z1}: */
            x = avg + dev * z;
            x1 = avg + dev * z1;
            v_pre1 = (which == 0 ? v_pre/dev : v_pre);
          }
        else
          { /* PDF is a Dirac pulse, CDF is a step function: */
            /* Pick {x,x1} different from {avg}: */
            x = avg + dmax * z;
            x1 = avg + dmax * z1;
            if (which == 0)
              { v_pre1 = (x != avg ? 0.0 : +INF); }
            else if (which == 1)
              { v_pre1 = (x < avg ? 0.0 : (x > avg ? 1.0 : NAN)); }
            else if (which == 2)
              { if ((x < avg) && (x1 > avg))
                  { v_pre1 = +1.0; }
                else if ((x > avg) && (x1 < avg)) 
                  { v_pre1 = -1.0; }
                else if ((x == avg) || (x1 == avg))
                  { v_pre1 = NAN; }
                else
                  { v_pre1 = 0.0; }
              }
          }
        double v_cmp = NAN;
        if (which == 0)
          { fname = "PDF";
            v_cmp = gauss_distr_PDF(x, avg, dev);
          }
        else if (which == 1)
          { fname = "CDF";
            v_cmp = gauss_distr_CDF(x, avg, dev);
          }
        else if (which == 2)
          { fname = "integral";
            v_cmp = gauss_distr_integral(x, x1, avg, dev);
          }
        test_check(fname, x, x1, avg, dev, v_cmp, v_pre1);
      }
  }

void test_check
  ( char *fname, 
    double x, 
    double x1, 
    double avg, 
    double dev, 
    double v_cmp, 
    double v_pre
  )
  {
    double ln_v_pre = (v_pre == 0.0 ? -INF : log(v_pre));
    double ln_v_cmp = (v_cmp == 0 ? -INF : log(fabs(v_cmp)));
    double errln;
    if (isnan(v_cmp) && isnan(v_pre))
      { errln = 0.0; }
    else if (v_cmp == v_pre)
      { errln = 0.0; }
    else if ((ln_v_cmp < -740.0) && (v_pre < -740.0))
      { errln = 0.0; }
    else
      { errln = ln_v_cmp - ln_v_pre; }
    fprintf(stderr, "%-20s\n", fname);
    fprintf(stderr, "  avg =       %24.16e  dev =       %24.16e\n", avg, dev);
    fprintf(stderr, "  x =         %24.16e", x);
    if (! isnan(x1)) { fprintf(stderr, "  x1 =        %24.16e", x1); }
    fprintf(stderr, "\n");
    fprintf(stderr, "  v_cmp =     %24.16e  v_pre =     %24.16e\n", v_cmp, v_pre);
    fprintf(stderr, "  ln(v_cmp) = %+24.16e  ln(v_pre) = %+24.16e\n", ln_v_cmp, ln_v_pre);
    fprintf(stderr, "  errln = %+24.16f\n", errln);
    /* We cannot require 16 digits precision because {ln_v_cmp} goes through {log(exp(*))}: */
    if (fabs(errln) > 1.0e-11) { fprintf(stderr, "   ** logs do not match\n"); }
    fprintf(stderr, "\n");
  }
