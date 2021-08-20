/* See rn_test_tools.h. */
/* Last edited on 2021-08-18 15:22:42 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>

#include <rn_test_tools.h>

#define NO NULL

void rn_do_check_eq(double x, double y, int32_t *i, int32_t *j, char *msg, rn_LOCPARMS)
  { if (x != y)
      { fprintf(stderr, " **"); 
        if (i != NULL) { fprintf(stderr, " [%d]", *i); }
        if (j != NULL) { fprintf(stderr, " [%d]", *j); }
        fprintf(stderr, " %20.16e != %20.16e\n", x, y);
        programerror(msg, file, line, func);
      }
  }

void rn_do_check_eps(double x, double y, double eps, int32_t *i, int32_t *j, char *msg, rn_LOCPARMS)
  { if (fabs(x - y) > eps)
      { fprintf(stderr, " **"); 
        if (i != NULL) { fprintf(stderr, " [%d]", *i); }
        if (j != NULL) { fprintf(stderr, " [%d]", *j); }
        double dif = fabs(x - y);
        fprintf(stderr, " | %20.16e - %20.16e | =  %20.16e > %20.16e\n", x, y, dif, eps);
        programerror(msg, file, line, func);
      }
  }

void rn_test_rot_axis(int32_t n, double *a, int32_t i, int32_t j, double ang, double *r, char *msg)
  { 
    assert((i >= 0) && (i < n));
    assert((j >= 0) && (j < n) && (i != j));
    double ca = cos(ang), sa = sin(ang);
    for (int32_t k = 0; k < n; k++)
      { double rrk;
        if (k == i) 
          { rrk = + ca*a[k] - sa*a[j]; }
        else if (k == j)
          { rrk = + sa*a[i] + ca*a[k]; }
        else
          { rrk = a[k]; }
        rn_check_eq(r[k], rrk, &i, &j, msg);
      }
  }
