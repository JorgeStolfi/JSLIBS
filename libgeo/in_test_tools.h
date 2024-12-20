/* in_test_tools.h --- tools for testing the {libgeo} functions. */
/* Last edited on 2024-12-05 10:27:42 by stolfi */
/* Created 2005-07-20 by J. Stolfi. */

#ifndef in_test_tools_H
#define in_test_tools_H

#include <stdio.h>
#include <stdint.h>

#include <affirm.h>

#define in_LOCPARMS  const char *file, int32_t line, const char *func
#define in_LOCARGS   __FILE__, __LINE__, __FUNCTION__

/* NUMERICAL CHECKS FOR INTEGER PROCEDURES
   
   Each procedure in this section is a variant of {affirm} that tests
   some condition of numerical values. If the test fails, the
   procedure prints the offending values, and aborts with an error
   message. In that case, if the {int32_t *} arguments {i} and/or {j} are
   not NULL, the procedure also prints the values of {*i} and/or
   {*j}. */

#define in_check_eq(x, y, i, j, msg)  in_do_check_eq(x, y, i, j, msg, in_LOCARGS)
  /* Checks that the {double}s {x} and {y} are exactly equal. */

void in_do_check_eq(int64_t x, int64_t y, uint32_t *i, uint32_t *j, char *msg, in_LOCPARMS);

#endif
