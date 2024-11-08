#define PROG_NAME "test_gauss_distr"
#define PROG_DESC "test of {gauss_distr.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-11-05 21:05:56 by stolfi */ 
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

void test_all(double x, double avg, double dev);
  /* Tests the {gauss_distr_PDF(x, avg, dev)},
    {gauss_distr_CDF(x, avg, dev)} and {gauss_distr_integral(x0, x1, avg, dev)}.
    . */


void test_gauss_PDF(double x, double avg, double dev);
  /* Tests {gauss_distr_PDF(x, avg, dev)} */
    
void test_gauss_CDF(double x, double avg, double dev);
  /* Tests {gauss_distr_CDF(x, avg, dev)}. */

void test_gauss_integral(double x0, double x1, double avg, double dev);
  /* Tests {gauss_distr_integral(x0, x1, avg, dev)} */

void test_check
  ( char *fname, 
    double x0, 
    double x1, 
    double avg,
    double dev,
    double v_cmp,
    double v_pre
  );
  /* Checks the computed result {v_cmp} of function
    {gauss_distr_PDF(x0,avg,dev)} or {gauss_distr_CDF(x0,avg,dev)} (or
    {gauss_distr_{fname}(x0,x1,avg,dev)} if {x1} is not {NAN}) with the
    expected value {v_pre}. */

double test_arg(int32_t k, int32_t N, double sgn, double vmin, double vmax);
  /* Returns a test value from a series.
    
    If {k} is in {-N .. +N}, returns an element of the geometric series
    {sgn*vmin*(vmax/vmin)^k/N}.  The values {vmin} and {vmax} must be finite an 
    positive.
    
    If {|k|} is {N+2} or {N+1}, returns {0}, or {sgn*INF},
    respectively. */

int32_t main(int32_t argn, char **argv);

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  {
    test_basics();
    return 0;
  }

double test_arg(int32_t k, int32_t N, double sgn, double vmin, double vmax)
  { bool_t debug = TRUE;
    if (debug) { fprintf(stderr, "    k = %d N = %d sgn = %+.0f vmin = %24.16e vmax = %24.16e\n", k, N, sgn, vmin, vmax); }
    assert((sgn == -1) || (sgn == +1));
    assert(isfinite(vmin) && vmin > 0);
    assert(isfinite(vmax) && vmax >= vmin);
    double v;
    if (abs(k) == N+2)
      { v = 0.0; }
    else if (abs(k) == N+1)
      { v = sgn*INF; }
    else if (vmax == vmin)
      { return sgn*vmin; }
    else
      { double fr = fabs(((double)k+N)/((double)2*N));
        v = sgn * vmin * pow(vmax/vmin, fr);
      }
    if (debug) { fprintf(stderr, "    v = %24.16e\n", v); }
    assert(! isnan(v));
    return v;
  }

void test_basics(void)
  {
    bool_t debug = TRUE;

    double avg = 20.0;
 
    /* Deviations to consider: */
    int32_t Nd = 3;
    double dmin = 1.0e-6, dmax = 1.0e+6;
   
    /* Test values for {u}: */
    int32_t Nz = 5;
    double umin = 1.0e-6, umax = 1.0e+6;

    for (int32_t kd = -Nd; kd <= Nd+2; kd++)
      { double dev = test_arg(kd, Nd, +1, dmin, dmax); 
        double sgn0 = +1;
        for (int32_t kz0 = -Nz; kz0 <= Nz+2; kz0++)
          { double u0 = test_arg(kz0, Nz, sgn0, umin, umax);
            double x0 = avg + (isfinite(dev) && (dev != 0) ? dev : 1) * u0;
            if (debug) { fprintf(stderr, "    u0 = %24.16e  x0 = %24.16e\n", u0, x0); }
            test_gauss_PDF(x0, avg, dev);
            test_gauss_CDF(x0, avg, dev);
            double sgn1 = +1;
            for (int32_t kz1 = -Nz; kz1 <= Nz+2; kz1++)
              { double u1 = test_arg(kz1, Nz, sgn1, umin, umax);
                double x1 = avg + (isfinite(dev) && (dev != 0) ? dev : 1) * u1;
                if (debug) { fprintf(stderr, "    u1 = %24.16e  x1 = %24.16e\n", u1, x1); }
                test_gauss_integral(x0, x1, avg, dev);
                sgn1 = -sgn1;
              }
            sgn0 = -sgn0;
          }
      }
   }
   
void test_gauss_PDF(double x, double avg, double dev)
  { fprintf(stderr, "    > --- %s x = %24.16e ---\n", __FUNCTION__, x);
    demand(! isnan(x), "bad {x0}");
    demand(isfinite(avg), "bad {avg}");
    demand(! isnan(dev), "bad {dev}");
    double v_pre; /* Predicted value for this {x1}. */
    if (dev == 0.0)
      { /* PDF is a Dirac pulse: */
        v_pre = (x != avg ? 0.0 : +INF);
      }
    else
      { /* PDF is bell-shaped: */
        double u = (fabs(x) == INF ? x : (x - avg)/dev);
        v_pre = exp(-u*u/2)/sqrt(2.0*M_PI)/dev;
      }
    double v_cmp = gauss_distr_PDF(x, avg, dev);
    test_check("PDF", x, NAN, avg, dev, v_cmp, v_pre);
    fprintf(stderr, "    < --- %s ---\n", __FUNCTION__);
  }
   
void test_gauss_CDF(double x, double avg, double dev)
  { fprintf(stderr, "    > --- %s x = %24.16e ---\n", __FUNCTION__, x);
    demand(! isnan(x), "bad {x0}");
    demand(isfinite(avg), "bad {avg}");
    demand(! isnan(dev), "bad {dev}");
    double v_pre; /* Predicted value for this {x1}. */
    if (dev == 0.0)
      { /* CDF is a step at {u==0}: */
        v_pre = (x < avg ? 0.0 : ( x > avg ? 1.0 : NAN));
      }
    else if (dev == +INF)
      { /* CDF is zero everywhere, but {NAN} at {+INF}: */
        v_pre = (x == +INF ? 1.0 : 0.0);
      }
    else
      { /* CDF is sigmoid: */
        double u = (fabs(x) == INF ? x : (x - avg)/dev);
        v_pre = 0.5*(1+ erf(u/M_SQRT2));;
      }
    double v_cmp = gauss_distr_CDF(x, avg, dev);
    test_check("CDF", x, NAN, avg, dev, v_cmp, v_pre);

    fprintf(stderr, "    < --- %s ---\n", __FUNCTION__);
  }
   
void test_gauss_integral(double x0, double x1, double avg, double dev)
  { fprintf(stderr, "    > --- %s x0 = %24.16e x1 = %24.16e ---\n", __FUNCTION__, x0, x1);

    bool_t debug = TRUE;

    demand(! isnan(x0), "bad {x0}");
    demand(! isnan(x1), "bad {x1}");
    demand(isfinite(avg), "bad {avg}");
    demand(! isnan(dev), "bad {dev}");
    double v_pre; /* Predicted value for this {x0,x1}. */
    if (dev == 0.0)
      { /* CDF is a step function at {avg}: */
        if ((x0 < avg) && (x1 > avg))
          { v_pre = +1.0; }
        else if ((x0 > avg) && (x1 < avg)) 
          { v_pre = -1.0; }
        else if ((x0 == avg) || (x1 == avg))
          { v_pre = NAN; }
        else
          { v_pre = 0.0; }
      }
    else if (dev == INF)
      { /* CDF is a zero function with unit integral: */
        if ((x0 == -INF) && (x1 == +INF))
          { v_pre = +1.0; }
        else if ((x0 == +INF) && (x1 == -INF))
          { v_pre = -1.0; }
        else if (x0 == x1)
          { v_pre = 0.0; }
        else if ((fabs(x0) == INF) || (fabs(x1) == INF))
          { v_pre = NAN; }
        else
          { v_pre = 0.0; }
      }
    else
      { /* CDF is a sigmoid: */
        double u0 = (fabs(x0) == INF ? x0 : (x0 - avg)/dev);
        double u1 = (fabs(x1) == INF ? x1 : (x1 - avg)/dev);
        if (fabs(u1 - u0) < 1.0e-4)
          { /* Numerical integration of PDF: */
            double g0 = gauss_distr_PDF(x0, avg, dev);
            double gm = gauss_distr_PDF(0.5*(x0+x1), avg, dev);
            double g1 = gauss_distr_PDF(x1, avg, dev);
            if (debug) { fprintf(stderr, "        g0 = %24.16e  gm = %24.16e  g1 = %24.16e  u1-u0 = %24.16e\n", g0, gm, g1, u1-u0); }
            v_pre = (0.25*g0 + gm + 0.25*g1)*(2.0/3.0)*(x1 - x0);
          }
        else
          { /* Define {ei = erf(ui) = ai + ci} where {ai} is {±1}. */
            double e0 = 0.5*erf(fabs(u0)/M_SQRT2), a0 = 0.5, c0 = -0.5*erfc(fabs(u0)/M_SQRT2); 
            if (u0 < 0) { e0 = -e0; a0 = -a0; c0 = -c0; }
            if (debug) { fprintf(stderr, "        u0 = %24.16e  e0 = %24.16e  c0 = %24.16e\n", u0, e0, c0); }
            double e1 = 0.5*erf(fabs(u1)/M_SQRT2), a1 = 0.5, c1 = -0.5*erfc(fabs(u1)/M_SQRT2);
            if (u1 < 0) { e1 = -e1; a1 = -a1; c1 = -c1; }
            if (debug) { fprintf(stderr, "        u1 = %24.16e  e1 = %24.16e  c1 = %24.16e\n", u1, e1, c1); }
            /* The integral is {e1-e0 = e1-(a0+c0) = (a1+c1)-e0 = (a1+c1)-(a0+c0)}: */
            if (fabs(e0) < fabs(c0)) 
              { if (fabs(e1) < fabs(c1))
                  { v_pre = e1-e0; }
                else
                  { v_pre = a1+(c1-e0); }
              }
            else
              { if (fabs(e1) < fabs(c1))
                  { v_pre = (e1-c0) - a0; }
                else
                  { v_pre = (c1-c0) + (a1-a0); }
              }
          }
      }
    double v_cmp = gauss_distr_integral(x0, x1, avg, dev);
    test_check("integral", x0, x1, avg, dev, v_cmp, v_pre);
    fprintf(stderr, "    < --- %s ---\n", __FUNCTION__);
  }

void test_check
  ( char *fname, 
    double x0, 
    double x1, 
    double avg, 
    double dev, 
    double v_cmp, 
    double v_pre
  )
  {
    fprintf(stderr, "      %-20s\n", fname);
    double u0, u1;
    if (dev == 0)
      { u0 = (isnan(x0) ? NAN : (x0 == avg ? 0.0 : (x0 < avg ? -INF : +INF)));
        u1 = (isnan(x1) ? NAN : (x1 == avg ? 0.0 : (x1 < avg ? -INF : +INF)));
      }
    else
      { u0 = (isnan(x0) ? NAN : (fabs(x0) == INF ? x0 : (x0 - avg)/dev));
        u1 = (isnan(x1) ? NAN : (fabs(x1) == INF ? x1 : (x1 - avg)/dev));
      }
    if (isnan(x1))
      { fprintf(stderr, "         x =          %24.16e  u =         %24.16e\n", x0, u0); }
    else
      { fprintf(stderr, "         x0 =         %24.16e  u0 =        %24.16e\n", x0, u0);
        fprintf(stderr, "         x1 =         %24.16e  u1 =        %24.16e\n", x1, u1);
      }
    fprintf(stderr, "        avg =         %24.16e  dev =       %24.16e\n", avg, dev);
    fprintf(stderr, "        v_cmp =       %24.16e  v_pre =     %24.16e\n", v_cmp, v_pre);
    
    bool_t ok = TRUE;
    if (isnan(v_pre) && isnan(v_cmp))
      { ok = TRUE; }
    else if (isnan(v_pre) || isnan(v_cmp))
      { fprintf(stderr, "         ** inconsistent {NAN}s\n");
        ok = FALSE;
      }
    else if (v_pre == v_cmp)
      { ok = TRUE; }
    else if ((v_pre == 0) || (v_cmp == 0))
      { fprintf(stderr, "         ** inconsistent zeros\n");
        ok = FALSE;
      }
    else if ((v_pre > 0) != (v_cmp > 0))
      { fprintf(stderr, "         ** different signs\n");
        ok = FALSE;
      }
    else
      { double ln_v_pre = (v_pre == 0.0 ? -INF : log(fabs(v_pre)));
        double ln_v_cmp = (v_cmp == 0 ? -INF : log(fabs(v_cmp)));
        double errln;
        if (v_cmp == v_pre)
          { errln = 0.0; }
        else if ((ln_v_cmp < -740.0) && (v_pre < -740.0))
          { errln = 0.0; }
        else
          { errln = ln_v_cmp - ln_v_pre; }
        fprintf(stderr, "        ln(|v_cmp|) = %+24.16e  ln(|v_pre|) = %+24.16e\n", ln_v_cmp, ln_v_pre);
        fprintf(stderr, "        errln =   %+24.16f\n", errln);
        /* We cannot require 16 digits precision because {ln_v_cmp} goes through {log(exp(*))}: */
        if (fabs(errln) > 5.0e-11) 
          { fprintf(stderr, "         ** logs do not match\n");
            ok = FALSE;
          }
      }
    fprintf(stderr, "\n");
    demand(ok, "failed");
  }
