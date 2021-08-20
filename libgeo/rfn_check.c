/* See rfn_check.h. */
/* Last edited on 2021-08-18 16:49:46 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <rn_check.h>

#include <rfn_check.h>

#define NO NULL

void rfn_check_rot_axis(int32_t n, float *a, int32_t i, int32_t j, double ang, float *r, char *msg)
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
        rfn_check_eq(r[k], rrk, &i, &j, msg);
      }
  }
