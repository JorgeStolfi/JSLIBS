/* test_rf3 --- test program for rf3.h, rf3x3.h  */
/* Last edited on 2024-11-20 21:30:43 by stolfi */

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <affirm.h>
#include <jsrandom.h>

#include <rf3.h>
#include <rfn_test_tools.h>

#define N 3
#define NO NULL

/* Internal prototypes */

int32_t main (int32_t argc, char **argv);
void test_rf3(bool_t verbose);

int32_t main (int32_t argc, char **argv)
  { srand(1993);
    srandom(1993);

    for (uint32_t i = 0;  i < 100; i++) { test_rf3(i < 3); }
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_rf3(bool_t verbose)
  { rf3_t a, b, d;
    double r, s, t;
    double rr, ss;

    if (verbose)
      { fprintf(stderr,
          "sizeof(rf3_t) = %lu  %d*sizeof(double) = %lu\n",
          sizeof(rf3_t), N, N*sizeof(double)
        );
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf3_throw_cube ---\n"); }
    a = rf3_throw_cube();
    for (uint32_t i = 0;  i < N; i++)
      { for (uint32_t j = 0;  j < i; j++)
          { affirm(a.c[i] != a.c[j], "rf3_throw_cube probable error(1)"); } 
        /* Check whether there are more than 8 nonzero bits: */
        double vv = a.c[i]*256.0;
        affirm(vv != floor(vv), "rf3_throw error(3)"); 
        affirm((a.c[i] > -1.0) && (a.c[i] < 1.0), "rf3_throw error(2)"); 
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf3_add ---\n"); }
    a = rf3_throw_cube();
    b = rf3_throw_cube();
    d = rf3_add(&a, &b);
    for (uint32_t i = 0;  i < N; i++)
      { rfn_test_tools_check_eq(d.c[i],a.c[i] + b.c[i], NO, NO, "rf3_add error"); }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf3_sub ---\n"); }
    a = rf3_throw_cube();
    b = rf3_throw_cube();
    d = rf3_sub(&a, &b);
    for (uint32_t i = 0;  i < N; i++)
      { rfn_test_tools_check_eq(d.c[i],a.c[i] - b.c[i], NO, NO, "rf3_sub error"); }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf3_scale ---\n"); }
    s = drandom();
    a = rf3_throw_cube();
    d = rf3_scale(s, &a);
    for (uint32_t i = 0;  i < N; i++)
      { float zi = (float)(s*a.c[i]);
        rfn_test_tools_check_eq(d.c[i],zi, NO, NO, "rf3_scale error(1)");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf3_mix ---\n"); }
    s = drandom();
    t = drandom();
    a = rf3_throw_cube();
    b = rf3_throw_cube();
    d = rf3_mix(s, &a, t, &b);
    for (uint32_t i = 0;  i < N; i++)
      { float ddi = (float)(s * a.c[i] + t * b.c[i]);
        rfn_test_tools_check_eq(d.c[i],ddi, NO, NO, "rf3_mix error");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf3_rot_axis ---\n"); }
    { a = rf3_throw_cube();
      uint32_t i = uint32_abrandom(0, N-1);
      uint32_t j = uint32_abrandom(0, N-2); if (j >= i) { j++; }
      double ang = 2.1*M_PI*drandom();
      d = rf3_rot_axis(&a, i, j, ang);
      rfn_test_tools_check_rot_axis(N, a.c, i, j, ang, d.c, "rf3_rot_axis error");
    }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf3_norm ---\n"); }
    a = rf3_throw_cube();
    r = rf3_norm(&a);
    ss = 0.0;
    for (uint32_t i = 0;  i < N; i++)
      { double ai = fabs(a.c[i]);
        ss += ai*ai; 
      }
    rr = sqrt(ss);
    rfn_test_tools_check_eps(r,rr,0.0000001 * rr, NO, NO, "rf3_norm error");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf3_max ---\n"); }
    if (verbose) { fprintf(stderr, "!! rf3_max NOT TESTED\n"); }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf3_rot_gen ---\n"); }
    if (verbose) { fprintf(stderr, "!! rf3_rot_gen NOT TESTED\n"); }
  }
