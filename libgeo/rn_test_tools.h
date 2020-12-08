/* rn_test_tools.h --- tools for testing the {libgeo} functions. */
/*
  Last edited on 2007-07-21 14:52:32 by stolfi
  Created 2005-07-20 by J. Stolfi.
*/

#ifndef rn_test_tools_H
#define rn_test_tools_H

#include <stdio.h>
#include <affirm.h>

#define rn_LOCPARMS  const char *file, int line, const char *func
#define rn_LOCARGS   __FILE__, __LINE__, __FUNCTION__

/* NUMERICAL CHECKS 
   
   Each procedure in this section is a variant of {affirm} that tests
   some condition of numerical values. If the test fails, the
   procedure prints the offending values, and aborts with an error
   message. In that case, if the {int *} arguments {i} and/or {j} are
   not NULL, the procedure also prints the values of {*i} and/or
   {*j}. */

#define rn_check_eq(x, y, i, j, msg)  rn_do_check_eq(x, y, i, j, msg, rn_LOCARGS)
  /* Checks that the {double}s {x} and {y} are exactly equal. */

#define rn_check_eps(x, y, eps, i, j, msg)  rn_do_check_eps(x, y, eps, i, j, msg, rn_LOCARGS)
  /* Checks that {double}s {x} and {y} satisfy {|x-y| <= eps}. */

void rn_do_check_eq(double x, double y, int *i, int *j, char *msg, rn_LOCPARMS);
void rn_do_check_eps(double x, double y, double eps, int *i, int *j, char *msg, rn_LOCPARMS);

/* OTHER TOOLS */

#endif
