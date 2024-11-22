/* See rfn_test_tools.h. */
/* Last edited on 2024-11-20 13:02:16 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <rn_test_tools.h>

#include <rfn_test_tools.h>

#define NO NULL

void rfn_test_tools_do_check_eq(double x, double y, uint32_t *i, uint32_t *j, char *msg, rfn_test_tools_LOCPARMS)
  { if (x != y)
      { fprintf(stderr, " **"); 
        if (i != NULL) { fprintf(stderr, " [%d]", *i); }
        if (j != NULL) { fprintf(stderr, " [%d]", *j); }
        fprintf(stderr, " %20.16e != %20.16e\n", x, y);
        programerror(msg, file, line, func);
      }
  }

void rfn_test_tools_do_check_eps(double x, double y, double eps, uint32_t *i, uint32_t *j, char *msg, rfn_test_tools_LOCPARMS)
  { double diff = fabs(x - y);
    if (diff > eps)
      { fprintf(stderr, " **"); 
        if (i != NULL) { fprintf(stderr, " [%d]", *i); }
        if (j != NULL) { fprintf(stderr, " [%d]", *j); }
        fprintf(stderr, " | %20.16e - %20.16e | =  %20.16e > %20.16e\n", x, y, diff, eps);
        programerror(msg, file, line, func);
      }
  }
  
void rfn_test_tools_check_rot_axis(uint32_t n, float *a, uint32_t i, uint32_t j, double ang, float *r, char *msg)
  { 
    assert((i >= 0) && (i < n));
    assert((j >= 0) && (j < n) && (i != j));
    double ca = cos(ang), sa = sin(ang);
    for (int32_t k = 0; k < n; k++)
      { float rrk;
        if (k == i) 
          { rrk = (float)(+ ca*a[k] - sa*a[j]); }
        else if (k == j)
          { rrk = (float)(+ sa*a[i] + ca*a[k]); }
        else
          { rrk = a[k]; }
        rfn_test_tools_check_eq(r[k], rrk, &i, &j, msg);
      }
  }
