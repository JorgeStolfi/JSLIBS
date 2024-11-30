/* test_rf2 --- test program for rf2.h, rf2x2.h  */
/* Last edited on 2024-11-20 13:29:34 by stolfi */

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <affirm.h>
#include <jsrandom.h>
#include <rfn_test_tools.h>

#include <rf2.h>

#define N 2
#define NO NULL

/* Internal prototypes */

int32_t main (int32_t argc, char **argv);
void test_rf2(bool_t verbose);

int32_t main (int32_t argc, char **argv)
  {
    srand(1993);
    srandom(1993);

    for (uint32_t i = 0;  i < 100; i++) { test_rf2(i < 3); }
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_rf2(bool_t verbose)
  {
    rf2_t a, d;
    double s;

    if (verbose)
      { fprintf(stderr,
          "sizeof(rf2_t) = %lu  %d*sizeof(double) = %lu\n",
          sizeof(rf2_t), N, N*sizeof(double)
        );
      }
    float ff = frandom();
    assert((ff >= 0.0) && (ff < 1.0));

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf2_throw_cube ---\n"); }
    a = rf2_throw_cube();
    for (uint32_t i = 0;  i < N; i++)
      { for (uint32_t j = 0;  j < i; j++)
          { affirm(a.c[i] != a.c[j], "rf2_throw_cube probable error(1)"); } 
        affirm((a.c[i] > -1.0) && (a.c[i] < 1.0), "rf2_throw_cube error(2)"); 
        /* Check whether there are more than 8 nonzero bits: */
        double vv = a.c[i]*256.0;
        affirm(vv != floor(vv), "rf2_throw_cube error(3)"); 
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf2_scale ---\n"); }
    s = drandom();
    a = rf2_throw_cube();
    d = rf2_scale(s, &a);
    for (uint32_t i = 0;  i < N; i++)
      { float zi = (float)(s*a.c[i]);
        rfn_test_tools_check_eq(d.c[i],zi, NO, NO, "rf2_scale error(1)");
      }
  }

