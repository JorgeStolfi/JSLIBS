/* rn_test_tools.h --- tools for testing the {libgeo} functions. */
/* Last edited on 2024-12-05 10:28:51 by stolfi */
/* Created 2005-07-20 by J. Stolfi. */

#ifndef rn_test_tools_H
#define rn_test_tools_H

#include <stdio.h>
#include <stdint.h>

#include <affirm.h>
#include <rn.h>

#define rn_test_tools_LOCPARMS  const char *file, int32_t line, const char *func
#define rn_test_tools_LOCARGS   __FILE__, __LINE__, __FUNCTION__

/* NUMERICAL CHECKS 
   
   Each procedure in this section is a variant of {affirm} that tests
   some condition of numerical values. If the test fails, the
   procedure prints the offending values, and aborts with an error
   message. In that case, if the {uint32_t *} arguments {i} and/or {j} are
   not NULL, the procedure also prints the values of {*i} and/or
   {*j}. */

#define rn_test_tools_check_eq(x, y, i, j, msg) \
  rn_test_tools_do_check_eq(x, y, i, j, msg, rn_test_tools_LOCARGS)
  /* Checks that the {double}s {x} and {y} are exactly equal. */

#define rn_test_tools_check_eps(x, y, eps, i, j, msg) \
  rn_test_tools_do_check_eps(x, y, eps, i, j, msg, rn_test_tools_LOCARGS)
  /* Checks that {double}s {x} and {y} satisfy {|x-y| <= eps}. */

void rn_test_tools_do_check_eq(double x, double y, uint32_t *i, uint32_t *j, char *msg, rn_test_tools_LOCPARMS);
void rn_test_tools_do_check_eps(double x, double y, double eps, uint32_t *i, uint32_t *j, char *msg, rn_test_tools_LOCPARMS);

/* OTHER TOOLS */

void rn_test_tools_check_all_different(uint32_t n, double *a, char *msg);
  /* Compares all pairs of elements of {a} and bombs out if
    two are equal.  Used to check vectors randomly sampled
    from a distribution where any such identity is extremely unlikely,
    and most probably indicates an error in the generator. */

void rn_test_tools_check_rot_axis(uint32_t n, double *a, uint32_t i, uint32_t j, double ang, double *r, char *msg);
  /* Checks whether {r} is indeed the result of {rn_rot_axis(n, a, i, j, ang, r)}. */

#endif
