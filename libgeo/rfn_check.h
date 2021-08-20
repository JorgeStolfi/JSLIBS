/* rfn_check.h --- tools for testing the {libgeo} functions (single-precision). */
/* Last edited on 2021-08-18 16:50:45 by stolfi */
/* Created 2005-07-20 by J. Stolfi. */

#ifndef rfn_check_H
#define rfn_check_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <affirm.h>

#include <rn_test_tools.h>
#include <rfn_check.h>

#define rfn_LOCPARMS  const char *file, int32_t line, const char *func
#define rfn_LOCARGS   __FILE__, __LINE__, __FUNCTION__

/* NUMERICAL CHECKS 
   
   Each procedure in this section is a variant of {affirm} that tests
   some condition of numerical values. If the test fails, the
   procedure prints the offending values, and aborts with an error
   message. In that case, if the {int32_t *} arguments {i} and/or {j} are
   not NULL, the procedure also prints the values of {*i} and/or
   {*j}. */

/* OTHER TOOLS */

void rfn_check_rot_axis(int32_t n, float *a, int32_t i, int32_t j, double ang, float *r, char *msg);
  /* Checks whether {r} is indeed the result of {rfn_rot_axis(n, a, i, j, ang, r)}. */

#endif
