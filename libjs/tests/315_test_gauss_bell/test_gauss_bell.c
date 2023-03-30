#define PROG_NAME "test_gauss_bell"
#define PROG_DESC "test of {gauss_bell.h}"
#define PROG_VERS "1.0"

/* Last edited on 2023-03-18 11:29:33 by stolfi */ 
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

#include <gauss_bell.h>

#include <bool.h>
#include <jsmath.h>
#include <jsfile.h>
#include <affirm.h>

void test_basics(void);
  /* Runs tests on the functions of the interface {gauss_bell.h}. */

void test_eval(double z, double ln_v_exp);
  /* Tests {gauss_bell_eval(avg+dev*z, avg, dev)}
    for various values of {avg} and {dev}. Compares the natural log of the
    result with {ln_v_exp}. */
    
void test_check(char *fname, double x, double avg, double dev, double v_cmp, double ln_v_exp);
  /* Compares the log of the computed result {v_cmp} of some function of argument {x}
    with {ln_v_exp}.  The string {fname} identifies
    which function was called. */

double test_arg(int32_t k, int32_t N, double vmin, double vmax);
  /* Returns a test value from a series. If {k} is in {-N .. +N},
    returns {(-1)^k*vmin*(vmax/vmin)^k/N}. If {k} is in {-N-3}, {-N-2}, or {-N-1},
    returns {-INF}, {0}, or {+INF}, respectively. */

int32_t main(int32_t argn, char **argv);

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  {
    test_basics();
    return 0;
  }

double test_arg(int32_t k, int32_t N, double vmin, double vmax)
  { double v;
    if (k < -N)
      { v = ((N + k) + 2.0) * INF; }
    else
      { double s = 1.0 - 2*(k % 2);
        double fr = ((double)k)/((double)N);
        v = s * vmin * pow(vmax/vmin, fr);
      }
    return v;
  }

void test_basics(void)
  {
    test_eval(-INF, -INF);
    test_eval(+INF, -INF);
    test_eval(0.0, 0.0);
    int32_t Nz = 10;
    double zmin = 1.0e-6, zmax = 1.0e+6;
    for (int32_t kz = -3; kz <= Nz; kz++)
      { double z = test_arg(kz, Nz, zmin, zmax);
        double ln_v_exp = -z*z/2;
        test_eval(z, ln_v_exp);
      }
   }
   
void test_eval(double z, double ln_v_exp)
  { double avg = 20.0; 
    double dmin = 0.01, dmax = 10.0; /* Range of deviations to try. */
    int32_t Nd = 3; /* Number of deviations to try, besides {±INF} and zero. */
    for (int32_t kd = -3; kd <= Nd; kd++)
      { double dev = test_arg(kd, Nd, dmin, dmax);
        double x; /* Argument of generic bell function. */
        double ln_vx_exp; /* Log of expected result of {gauss_bell_eval(x, avg, dev)}. */
        if (dev == 0)
          { /* Pick a {x} argument different from {avg}: */
            x = avg + dmax * z;
            ln_vx_exp = (z != 0.0 ? -INF : 0.0);
          }
        else if (fabs(dev) == +INF)
          { /* Pick a {x} argument different from {avg}: */
            x = avg + dmax * z;
            ln_vx_exp = 0.0;
          }
        else
          { /* Pick the {x} argument that corresponds to the standard {z}: */
            x = avg + dev * z;
            ln_vx_exp = ln_v_exp;
          }
        double vx_cmp = gauss_bell_eval(x, avg, dev);
        test_check("generic", x, avg, dev, vx_cmp, ln_vx_exp);
      }
  }

void test_check(char *fname, double x, double avg, double dev, double v_cmp, double ln_v_exp)
  {
    double ln_v_cmp = (v_cmp == 0 ? -INF : log(v_cmp));
    double v_exp = (ln_v_exp == -INF ? 0.0 : exp(ln_v_exp));
    double errln;
    if (isnan(ln_v_cmp) && isnan(ln_v_exp))
      { errln = 0.0; }
    else if (ln_v_cmp == ln_v_exp)
      { errln = 0.0; }
    else if ((ln_v_cmp < -740.0) && (ln_v_exp < -740.0))
      { errln = 0.0; }
    else
      { errln = ln_v_cmp - ln_v_exp; }
    fprintf(stderr, "%-20s\n", fname);
    fprintf(stderr, "  avg =       %24.16e  dev =       %24.16e\n", avg, dev);
    fprintf(stderr, "  x =         %24.16e\n", x);
    fprintf(stderr, "  v_cmp =     %24.16e  v_exp =     %24.16e\n", v_cmp, v_exp);
    fprintf(stderr, "  ln(v_cmp) = %+24.16e  ln(v_exp) = %+24.16e\n", ln_v_cmp, ln_v_exp);
    fprintf(stderr, "  errln = %+24.16f\n", errln);
    /* We cannot require 16 digits precision because {ln_v_cmp} goes through {log(exp(*))}: */
    if (isnan(errln) || (fabs(errln) > 1.0e-9)) { fprintf(stderr, "   ** some error\n"); }
    fprintf(stderr, "\n");
  }
