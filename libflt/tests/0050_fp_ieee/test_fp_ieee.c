/* A short program to check several detaisl of 
   IEEE floating-point arithmetic:
     - whether (-oo) + (+oo) = NaN in all rounding modes.
     - how negative zero behaves with operations.
*/
/* Last edited on 2016-12-26 18:43:59 by stolfilocal */

#define _GNU_SOURCE
#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <values.h>
#include <math.h>
#include <stdio.h>
 
// #if (defined(SunOS5))
// 
// #include <ieeefp.h>
// 
// #define ROUND_NEAR { fpsetround(FP_RN); }
// #define ROUND_UP   { fpsetround(FP_RP); }
// #define ROUND_DOWN { fpsetround(FP_RM); }
// 
// #endif

#include <fenv.h>
#include <fpu_control.h>

#define ROUND_NEAR { fesetround(FE_TONEAREST); }
#define ROUND_UP   { fesetround(FE_UPWARD); }
#define ROUND_DOWN { fesetround(FE_DOWNWARD); }

void test_func_rounding(char *name, double func(double x), double xmin, double xmax);
  /* Tests rounding of {func} for {x} between {xmin} and {xmax} in geometric progression. */

void prt_1_bool_t(char *title, int m, int p);
void prt_1_double(char *title, double m, double p);
void prt_2_bool_t(char *title, int mm, int mp, int pm, int pp);
void prt_2_double(char *title, double mm, double mp, double pm, double pp);

int main(int argc, char **argv)
{
  { /* Testing special values: */
    fprintf(stderr, "MAXDOUBLE = %24.16e\n", MAXDOUBLE);
    fprintf(stderr, "DBL_MAX =   %24.16e\n", DBL_MAX);
    fprintf(stderr, "HUGE_VAL =  %24.16e\n", HUGE_VAL);
    fprintf(stderr, "INF =       %24.16e\n", INFINITY);
    fprintf(stderr, "NAN =       %24.16e\n", NAN);
  }
  
  test_func_rounding("exp",  exp,   0.5, 0.5e003);
  test_func_rounding("sin",  sin,   0.5, 0.5e003);
  test_func_rounding("sqrt", sqrt,  0.5, 0.5e200);
  
  { /* Testing whether (-oo) + (+oo) = NaN in all rounding modes: */
    double one, zero, pinf, minf, x, y, z;
    one = 1.0;
    zero = 0.0;
    ROUND_UP;
    pinf = one/zero;
    minf = -pinf;
    ROUND_DOWN;
    x = minf + pinf;
    ROUND_UP;
    y = minf + pinf;
    ROUND_NEAR;
    z = minf + pinf;
    fprintf(stderr, "(-oo)+(+oo): DN = %g  NE = %g UP = %g\n", x, y, z);
  }
  
  { /* Testing the behavior of negative zero in various ops: */
    double mz, pz;
    int i;
    pz = 0.0;
    mz = -pz;
    for (i = 0; i < 3; i++)
      {
        if      (i == 0) { ROUND_NEAR; fprintf(stderr, "\n=== ROUND_NEAR ===\n\n"); }
        else if (i == 1) { ROUND_UP;   fprintf(stderr, "\n=== ROUND_UP ===\n\n"); }
        else if (i == 2) { ROUND_DOWN; fprintf(stderr, "\n=== ROUND_DOWN ===\n\n"); }
        
        prt_1_double("neg",  -mz, -pz);
        prt_1_double("sqrt", sqrt(mz), sqrt(pz));
        prt_1_double("ln",   log(mz), log(pz));
        prt_1_double("inv",  1.0/mz, 1.0/pz);
        prt_2_double("add",  mz+mz, mz+pz, pz+mz, pz+pz);
        prt_2_double("sub",  mz-mz, mz-pz, pz-mz, pz-pz);
        prt_2_double("mul",  mz*mz, mz*pz, pz*mz, pz*pz);

        prt_2_double("eql",  (mz==mz), (mz==pz), (pz==mz), (pz==pz));
        prt_2_double("gtr",  (mz>mz), (mz>pz), (pz>mz), (pz>pz));

      }        
  }
  return 0;
}

void test_func_rounding(char *name, double func(double x), double xmin, double xmax)
{ /* Testing rounding of {func}: */
  int N = 20;
  int bugs = 0;
  int i;
  for (i = 0; i <= N; i++)
    { ROUND_NEAR;
      double x = xmin * exp(i*(log(xmax) - log(xmin))/N);
      fprintf(stderr, "%s(%20.12e):\n", name, x);
      ROUND_NEAR;
      double ne = func(x);
      fprintf(stderr, "  NE = %20.12e\n", ne); 
      ROUND_DOWN;
      double dn = func(x);
      fprintf(stderr, "  DN = %20.12e (%20.12e)\n", dn, dn-ne); 
      ROUND_UP;
      double up = func(x);
      fprintf(stderr, "  UP = %20.12e (%20.12e)\n", up, up-ne);
      if ( ne > up ) { fprintf(stderr, "  ** NE > UP !\n"); bugs++; }
      if ( ne < dn ) { fprintf(stderr, "  ** NE < DN !\n"); bugs++; }
    }
  if (bugs > 0) { fprintf(stderr, "  ** %s does not round properly !\n", name); }
  fprintf(stderr, "\n");
}


void prt_1_bool_t(char *title, int m, int p)
{
  fprintf(stderr, "\n");
  fprintf(stderr, "%s\n", title);
  fprintf(stderr, "+----+----+----+\n");
  fprintf(stderr, "|    | -0 | +0 |\n");
  fprintf(stderr, "+----+----+----+\n");
  fprintf(stderr, "| -0 | %2d | %2d |\n", m, p);
  fprintf(stderr, "+----+----+----+\n");
}

void prt_1_double(char *title, double m, double p)
{
  fprintf(stderr, "\n");
  fprintf(stderr, "%s\n", title);
  fprintf(stderr, "+----+------+------+\n");
  fprintf(stderr, "|    |  -0  |  +0  |\n");
  fprintf(stderr, "+----+------+------+\n");
  fprintf(stderr, "| -0 | %4.0f | %4.0f |\n", m, p);
  fprintf(stderr, "+----+------+------+\n");
}

void prt_2_bool_t(char *title, int mm, int mp, int pm, int pp)
{
  fprintf(stderr, "\n");
  fprintf(stderr, "%s\n", title);
  fprintf(stderr, "+----+----+----+\n");
  fprintf(stderr, "|    | -0 | +0 |\n");
  fprintf(stderr, "+----+----+----+\n");
  fprintf(stderr, "| -0 | %2d | %2d |\n", mm, mp);
  fprintf(stderr, "+----+----+----+\n");
  fprintf(stderr, "| +0 | %2d | %2d |\n", pm, pp);
  fprintf(stderr, "+----+----+----+\n");
}

void prt_2_double(char *title, double mm, double mp, double pm, double pp)
{
  fprintf(stderr, "\n");
  fprintf(stderr, "%s\n", title);
  fprintf(stderr, "+----+------+------+\n");
  fprintf(stderr, "|    |  -0  |  +0  |\n");
  fprintf(stderr, "+----+------+------+\n");
  fprintf(stderr, "| -0 | %4.0f | %4.0f |\n", mm, mp);
  fprintf(stderr, "+----+------+------+\n");
  fprintf(stderr, "| +0 | %4.0f | %4.0f |\n", pm, pp);
  fprintf(stderr, "+----+------+------+\n");
}

